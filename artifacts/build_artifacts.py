"""Recompute v4 tables and figures from frozen JSONL records. No benchmarks run."""
from pathlib import Path
from collections import defaultdict, Counter
import csv
import hashlib
import json
import math
import os
import statistics as st

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
os.environ['MPLCONFIGDIR'] = str(ROOT / '_build' / 'matplotlib')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def read_json(path):
    return json.loads(path.read_text(encoding='utf-8'))

def read_rows(path):
    return [json.loads(line) for line in path.read_text(encoding='utf-8').splitlines() if line.strip()]

def write_json(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, ensure_ascii=False) + '\n', encoding='utf-8')

def write_csv(path, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('w', encoding='utf-8', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=list(dict.fromkeys(k for r in rows for k in r)))
        writer.writeheader()
        writer.writerows(rows)

def gm(values):
    return math.exp(st.mean(math.log(v) for v in values))

def check_close(actual, expected, context):
    assert math.isclose(actual, expected, rel_tol=1e-11, abs_tol=1e-12), (context, actual, expected)

def paired(rows, methods, keys, same_diagonal=False):
    groups = defaultdict(dict)
    for r in rows:
        assert r['status'] == 'ok', r
        key = tuple(r[k] for k in keys)
        assert r['method'] not in groups[key], key
        groups[key][r['method']] = r
    for key, by_method in groups.items():
        assert set(by_method) == set(methods), key
        vals = list(by_method.values())
        for field in ['graph_hash', 'answer'] + (['c', 'unique', 'shortcut'] if same_diagonal else []):
            assert len({r[field] for r in vals}) == 1, (key, field)
    return groups

# Verify frozen input hashes before using any data.
manifest = read_json(HERE / 'inputs_manifest.json')
for entry in manifest:
    assert hashlib.sha256((ROOT / entry['path']).read_bytes()).hexdigest() == entry['sha256'], entry

poly = HERE / 'theory' / 'data'
formal = sum((read_rows(poly / f'{s}.jsonl') for s in ['size', 'breadth', 'lower']), [])
assert len(formal) == 98
theory_groups = paired(formal, ['bh-ce', 'kw-ce'], ['suite', 'n', 'graph_seed', 'kind'])
theory_pairs = []
for key, rows in sorted(theory_groups.items()):
    bh, kw = rows['bh-ce'], rows['kw-ce']
    assert all(r['threads'] == 1 and r['full_determinization'] for r in rows.values())
    assert kw['layout'] == 'half'
    if bh['c'] == kw['c']:
        assert bh['unique'] == kw['unique']
    p = dict(suite='size' if key[0] == 'lower' else key[0], n=key[1], graph_seed=key[2], kind=key[3],
             bh_s=bh['total_s'], our_s=kw['total_s'], bh_over_our=bh['total_s']/kw['total_s'],
             bh_peak_mib=bh['rss_kib']/1024, our_peak_mib=kw['rss_kib']/1024,
             bh_choose_s=bh['choose_s'], our_setup_s=kw['setup_s'], our_choose_s=kw['choose_s'],
             bh_enumerate_s=bh['enumerate_s'], our_enumerate_s=kw['enumerate_s'],
             local_cache_bytes=kw['local_cache_bytes'], source_suite=key[0])
    theory_pairs.append(p)

small = read_rows(poly / 'dp_compare.jsonl')
assert len(small) == 189
small_groups = paired(small, ['bh-ce', 'kw-ce', 'dp-bit'], ['n', 'graph_seed', 'repetition'])
repeated = defaultdict(list)
reference = {(r['method'], r['n'], r['graph_seed']): r for r in formal if r['kind'] == 'dense'}
for r in small:
    assert r['threads'] == 1
    ref = reference['bh-ce', r['n'], r['graph_seed']]
    assert (r['answer'], r['graph_hash']) == (ref['answer'], ref['graph_hash'])
    if r['method'] == 'dp-bit':
        assert r['table_bytes'] == 8 * 2**(r['n']-1)
    else:
        ref = reference[r['method'], r['n'], r['graph_seed']]
        for field in ['c', 'unique', 'visits', 'centers']:
            assert r[field] == ref[field]
    repeated[r['n'], r['graph_seed'], r['method']].append(r)
small_pairs = []
for n, seed in sorted({(r['n'], r['graph_seed']) for r in small}):
    p = dict(n=n, graph_seed=seed, repetitions=3)
    for label, method in [('dp', 'dp-bit'), ('bh', 'bh-ce'), ('our', 'kw-ce')]:
        rows = repeated[n, seed, method]
        assert len(rows) == 3 and {r['repetition'] for r in rows} == {0, 1, 2}
        times = [r['total_s'] for r in rows]
        p.update({label+'_s': st.median(times), label+'_repeat_min_s': min(times),
                  label+'_repeat_max_s': max(times), label+'_peak_mib': max(r['rss_kib']/1024 for r in rows)})
    p.update(dp_over_our=p['dp_s']/p['our_s'], bh_over_our=p['bh_s']/p['our_s'],
             dp_table_mib=8*2**(n-1)/1024**2)
    small_pairs.append(p)

def summarize(rows, suite, methods):
    out = dict(suite=suite, n=rows[0]['n'], kind=rows[0].get('kind', 'dense'), graphs=len(rows))
    for method in methods:
        ts = [r[method+'_s'] for r in rows]
        out.update({method+'_s': st.median(ts), method+'_min_s': min(ts), method+'_max_s': max(ts),
                    method+'_peak_mib': max(r[method+'_peak_mib'] for r in rows)})
    out.update(bh_over_our=gm(r['bh_over_our'] for r in rows),
               our_wins_bh=sum(r['bh_over_our'] > 1 for r in rows))
    if 'dp' in methods:
        out.update(dp_over_our=gm(r['dp_over_our'] for r in rows),
                   our_wins_dp=sum(r['dp_over_our'] > 1 for r in rows),
                   dp_table_mib=rows[0]['dp_table_mib'])
    return out

theory_summary = []
for n in range(20, 45, 2):
    if n <= 32:
        rows = [r for r in small_pairs if r['n'] == n]
        r = summarize(rows, 'size', ['dp', 'bh', 'our'])
        r.update(repetitions=3, source='dp_compare.jsonl')
    else:
        rows = [r for r in theory_pairs if r['suite'] == 'size' and r['n'] == n]
        r = summarize(rows, 'size', ['bh', 'our'])
        r.update(repetitions=1, source='size.jsonl')
    theory_summary.append(r)
breadth = []
for n, kind in sorted({(r['n'], r['kind']) for r in theory_pairs if r['suite'] == 'breadth'}):
    rows = [r for r in theory_pairs if r['suite'] == 'breadth' and r['n'] == n and r['kind'] == kind]
    breadth.append(summarize(rows, 'breadth', ['bh', 'our']))

# Compare recomputed values against the archived publication summaries.
small_ref = {r['n']: r for r in read_json(poly/'dp_compare_summary.json')}
formal_ref = {(r['suite'], r['n'], r['kind']): r for r in read_json(poly/'summary.json')}
for r in theory_summary + breadth:
    if r['suite'] == 'size' and r['n'] <= 32:
        ref = small_ref[r['n']]
        for dst, src in [('dp_s','dp_median_s'),('bh_s','bh_median_s'),('our_s','kw_median_s'),
                         ('dp_over_our','dp_over_kw_geomean'),('bh_over_our','bh_over_kw_geomean')]:
            check_close(r[dst],ref[src],(r['n'],dst))
    else:
        ref = formal_ref[r['suite'],r['n'],r['kind']]
        for dst, src in [('bh_s','bh_s_median'),('our_s','kw_s_median'),('bh_over_our','speedup_geomean')]:
            check_close(r[dst],ref[src],(r['n'],dst))

pract = HERE/'practical32'/'data'
practical_records = read_rows(pract/'corrected.jsonl')
assert len(practical_records) == 276
practical_groups = paired(practical_records, ['bh-practical','kw-practical'],
                          ['suite','n','graph_seed','kind','algorithm_seed','repetition'], True)
practical_pairs = []
for key, rows in sorted(practical_groups.items()):
    bh, kw = rows['bh-practical'], rows['kw-practical']
    assert bh['threads'] == kw['threads'] == 32
    p = dict(suite=key[0],n=key[1],graph_seed=key[2],kind=key[3],algorithm_seed=key[4],repetition=key[5],
             bh_s=bh['total_s'],our_s=kw['total_s'],bh_over_our=bh['total_s']/kw['total_s'],
             bh_peak_mib=bh['rss_kib']/1024,our_peak_mib=kw['rss_kib']/1024,
             shortcut=bh['shortcut'],bh_kernel=bh['kernel'],source_mode=bh['source_mode'],
             bh_started_unix=bh['started_unix'],our_started_unix=kw['started_unix'])
    practical_pairs.append(p)
practical_summary=[]
for suite,n,kind in sorted({(r['suite'],r['n'],r['kind']) for r in practical_pairs}):
    all_rows=[r for r in practical_pairs if (r['suite'],r['n'],r['kind'])==(suite,n,kind)]
    active=[r for r in all_rows if not r['shortcut']]
    out=summarize(active or all_rows,suite,['bh','our'])
    out.update(pairs=len(all_rows),timed_pairs=len(active),shortcuts=len(all_rows)-len(active))
    if not active:
        out.update(bh_over_our=None,our_wins_bh=0)
    # Preserve group peak memory, also for shortcut rows.
    out['our_peak_gib']=max(r['our_peak_mib'] for r in all_rows)/1024
    practical_summary.append(out)
pr_ref={(r['suite'],r['n'],r['kind']):r for r in read_json(pract/'summary.json')}
for r in practical_summary:
    ref=pr_ref[r['suite'],r['n'],r['kind']]
    for dst,src in [('bh_s','bh_median_s'),('our_s','kw_median_s'),('our_peak_gib','kw_peak_gib')]:
        check_close(r[dst],ref[src],(r['suite'],r['n'],r['kind'],dst))
    if r['timed_pairs']:
        check_close(r['bh_over_our'],ref['speedup_geomean'],(r['suite'],r['n'],r['kind']))
    assert r['our_wins_bh']==ref['kw_wins']

for folder,names in [('theory',{'summary':theory_summary,'breadth':breadth,'pairs_original':theory_pairs,'pairs_dp':small_pairs}),
                     ('practical32',{'summary':practical_summary,'pairs':practical_pairs})]:
    for name,rows in names.items():
        write_csv(HERE/folder/f'{name}.csv',rows)
        write_json(HERE/folder/f'{name}.json',rows)

def tex_table(path, head, align, rows):
    body=['% Generated by artifacts/build_artifacts.py; times are seconds.',
          r'\begin{tabular}{'+align+'}',r'\toprule',head+r'\\',r'\midrule']
    body += [' & '.join(row)+r'\\' for row in rows]
    body += [r'\bottomrule',r'\end{tabular}']
    path.write_text('\n'.join(body)+'\n',encoding='utf-8')

tex_table(HERE/'theory'/'size_table.tex',
          r'$n$ & Graphs & Reps. & BH & Ours & Packed DP & BH/ours & DP/ours',
          'rrrrrrrr',[[str(r['n']),str(r['graphs']),str(r['repetitions']),f"{r['bh_s']:.5f}" if r['n']<28 else f"{r['bh_s']:.3f}",
                      f"{r['our_s']:.5f}" if r['n']<28 else f"{r['our_s']:.3f}",
                      f"{r['dp_s']:.5f}" if 'dp_s' in r else '--',f"{r['bh_over_our']:.3f}",
                      f"{r['dp_over_our']:.3f}" if 'dp_s' in r else '--'] for r in theory_summary])
sizes=[r for r in practical_summary if r['suite']=='size']
tex_table(HERE/'practical32'/'size_table.tex',
          r'$n$ & Graphs & BH32 & Ours32 & BH/ours & Wins & Peak GiB',
          'rrrrrrr',[[str(r['n']),str(r['pairs']),f"{r['bh_s']:.4f}",f"{r['our_s']:.4f}",
                     f"{r['bh_over_our']:.3f}",f"{r['our_wins_bh']}/{r['pairs']}",f"{r['our_peak_gib']:.3f}"] for r in sizes])
stage_rows=[]
for n in [36,40,42,44]:
    rows=[r for r in theory_pairs if r['suite']=='size' and r['n']==n]
    stage_rows.append([str(n)]+[f'{st.median(r[k] for r in rows):.3f}' for k in
                      ['bh_choose_s','our_setup_s','our_choose_s','bh_enumerate_s','our_enumerate_s']])
tex_table(HERE/'theory'/'stages_table.tex',
          r'$n$ & BH diagonal & Our cover & Our diagonal & BH enumeration & Our enumeration',
          'rrrrrr',stage_rows)

plt.rcParams.update({'font.family':'DejaVu Sans','font.size':11,'axes.spines.top':False,
                     'axes.spines.right':False,'pdf.fonttype':42,'ps.fonttype':42,'svg.fonttype':'none'})
colors={'bh':'#C0642D','our':'#196991','dp':'#57565B'}
labels={'bh':'BH13','our':'Ours (KW cover)','dp':'Packed DP (exponential space)'}
markers={'bh':'s','our':'o','dp':'^'}
def time_curves(ax,rows,methods):
    for method in methods:
        active=[r for r in rows if method+'_s' in r]
        x=[r['n'] for r in active]; y=[r[method+'_s'] for r in active]
        err=[[r[method+'_s']-r[method+'_min_s'] for r in active],
             [r[method+'_max_s']-r[method+'_s'] for r in active]]
        ax.errorbar(x,y,yerr=err,fmt='-'+markers[method],markersize=4,capsize=2,
                    color=colors[method],label=labels[method],linewidth=1.4)
    ax.set_yscale('log'); ax.set_xlabel('Vertices n');ax.set_ylabel('Complete calculation (s)')
    ax.grid(axis='y',alpha=.18);ax.legend(frameon=False,fontsize=10)
def save_fig(fig,name):
    fig.tight_layout(pad=1.4)
    destination=HERE/('theory' if name.startswith('theory') else 'practical32')/'figures'
    destination.mkdir(parents=True,exist_ok=True)
    for ext in ['pdf','svg','png']:
        fig.savefig(destination/f'{name}.{ext}',dpi=220,bbox_inches='tight')
    plt.close(fig)

fig,axs=plt.subplots(1,2,figsize=(8.4,3.5))
time_curves(axs[0],theory_summary,['bh','our','dp'])
axs[0].set_title('(a) One thread, deterministic');axs[0].set_xticks([20,24,28,32,36,40,44])
for method in ['bh','our','dp']:
    rows=[r for r in theory_summary if method+'_peak_mib' in r]
    axs[1].plot([r['n'] for r in rows],[r[method+'_peak_mib'] for r in rows],
                '-'+markers[method],color=colors[method],markersize=4,label=labels[method])
axs[1].set_yscale('log');axs[1].set_xlabel('Vertices n');axs[1].set_ylabel('Maximum process RSS (MiB)')
axs[1].set_title('(b) Measured memory');axs[1].grid(axis='y',alpha=.18)
axs[1].set_xticks([20,24,28,32,36,40,44]);axs[1].legend(frameon=False,fontsize=10)
save_fig(fig,'theory_comparison')
fig,axs=plt.subplots(1,2,figsize=(8.4,3.5))
time_curves(axs[0],sizes,['bh','our']);axs[0].set_title('(a) Both 32 threads, Las Vegas')
for method in ['bh','our']:
    axs[1].plot([r['n'] for r in sizes],[r[method+'_peak_mib']/1024 for r in sizes],
                '-'+markers[method],color=colors[method],markersize=4,label=labels[method])
axs[1].set_yscale('log');axs[1].set_xlabel('Vertices n');axs[1].set_ylabel('Maximum process RSS (GiB)')
axs[1].set_title('(b) Memory cost');axs[1].legend(frameon=False,fontsize=10);axs[1].grid(axis='y',alpha=.18)
for ax in axs:ax.set_xticks([32,36,40,44,48,50])
save_fig(fig,'practical32_comparison')

def markdown_table(headers, rows):
    return '\n'.join(['| '+' | '.join(headers)+' |','| '+' | '.join(['---']*len(headers))+' |']+
                     ['| '+' | '.join(map(str,row))+' |' for row in rows])
theory_md='''# 理论配置：单线程、完整确定化、多项式空间

Ours 表示本项目的 KW 覆盖算法。BH 与 Ours 每次执行完整条件期望，给定输入后没有随机选择。
压位 DP 是单线程、确定性的**指数空间外部参照**，不归为多项式空间。
所有时间与内存均来自既有远程 i9-14900K 数据；v4 只复算统计与绘图，没有重新计时。

20–32 点：每图每方法三次，先按图取中位数，再在图间取中位数；三种方法来自同轮重测。
34–44 点：原规模组每图一次，44 点仅一图。比值为逐图时间比的几何平均，不是中位时间之比。
DP 在 32 点以上没有本组实测值，空格不表示超时或预测。图中误差棒为输入间最小值至最大值，不是置信区间。

'''
theory_md+=markdown_table(['n','图数×重复','BH 秒','Ours 秒','DP 秒','BH/Ours','DP/Ours','DP 表 MiB'],
 [[r['n'],f"{r['graphs']}×{r['repetitions']}",f"{r['bh_s']:.5f}",f"{r['our_s']:.5f}",
 f"{r['dp_s']:.5f}" if 'dp_s' in r else '—',f"{r['bh_over_our']:.3f}",
 f"{r['dp_over_our']:.3f}" if 'dp_s' in r else '—',f"{r['dp_table_mib']:.0f}" if 'dp_s' in r else '—'] for r in theory_summary])
theory_md+='''

20 点由 DP 最快，22–30 点由 BH 最快，32 点 BH 与本项目接近（BH/Ours 为 1.023）。
本项目在 24–32 点的 15 张输入上全部快于 DP；32 点约快 10.70 倍，DP 主表为 16 GiB。
34–44 点的 16 张稠密图，以及 36/40 点其他密度和 tournament 的 12 张图，本项目全部快于 BH。
44 点 561.159 秒对 1636.349 秒，2.916 倍仅属于该单例。

时间优势主要来自自环确定化。44 点 BH 选自环 1353.023 秒，本项目 317.473 秒，另有构图 18.145 秒。
严格流式算法使用 O(n²) 位；计时实现加入受 n² 条目上限控制的局部缓存，保守为 O(n³ log n) 位。
约 9 MiB RSS 包含运行库开销，不是该渐近空间界的证明。

原两方法正式配对 49 组/98 次，另有训练 18 次；三方重测为已有 21 张图上的 189 次。
现有验证记录包含 4341 张小图，DP 对 20–32 点提供独立答案核验。34–44 点仅比较共享 P3 的两方法答案。
原始记录在 data/，可编译源码在 code/，协议、环境与原验证记录在 provenance/。
'''
theory_md+='\n'+markdown_table(['n','图族','图数','BH 秒','Ours 秒','BH/Ours','Ours 胜数'],
 [[r['n'],r['kind'],r['graphs'],f"{r['bh_s']:.4f}",f"{r['our_s']:.4f}",f"{r['bh_over_our']:.3f}",r['our_wins_bh']] for r in breadth])+'\n'
(HERE/'theory'/'RESULTS.md').write_text(theory_md,encoding='utf-8')
practical_md='''# 实际配置：双方 32 线程、Las Vegas、允许指数空间

Ours 表示本项目的随机 KW 覆盖、随机自环、分桶排序去重及 P3。
BH 用随机自环与一般图 Algorithm C，在二部图上使用专用 Algorithm B。
两者都返回精确奇偶答案。Ours 保存访问列表，使用指数空间；BH 仍使用多项式空间。

结果全部来自既有远程 14900K（8 P 核、16 E 核、32 硬件线程、约 125 GiB 内存）。
v4 没有重跑性能实验。138 个比较实例含 120 张不同图；133 组为补测 BH32 配已有 Ours32，
另 5 组复用已有完整 32/32 配对。两边在不同测量时段取得的记录不能称为全组交替重测。
CPU 频率与后台服务未锁定；来源规则提前固定，没有挑选较快秒数。

速度比为逐实例 BH/Ours 的几何平均；预处理直接拒绝的 9 组单列，不纳入核心速度比。
组合数据是 276 条算法记录，不是 276 次新增运行；其中只有 133 次新增 BH32。
全部成功，原始记录核对图哈希、自环、P2 计数和答案，没有不一致。

'''
for suite in ['size','density','structure','randomness']:
    practical_md+='\n## '+suite+'\n\n'
    rows=[r for r in practical_summary if r['suite']==suite]
    practical_md+=markdown_table(['n','图族','配对','捷径','BH 秒','Ours 秒','BH/Ours','Ours 胜数','Ours 峰值 GiB'],
     [[r['n'],r['kind'],r['pairs'],r['shortcuts'],f"{r['bh_s']:.5f}",f"{r['our_s']:.5f}",
       f"{r['bh_over_our']:.3f}" if r['timed_pairs'] else '—',r['our_wins_bh'],f"{r['our_peak_gib']:.3f}"] for r in rows])+'\n'
practical_md+='''
32–38 点不支持稳定优势，BH 在各规模多数输入和中位时间上更快。
40–50 点稠密规模组 24 张图全部由本项目更快，50 点三图的几何平均为 2.485 倍，
中位时间 154.469 秒对 380.831 秒，但本项目最大峰值约 38.165 GiB。
六张二部图全部由 BH 专用 Algorithm B 更快，40/44 点分别约快 5.36/6.60 倍。
因此不能声称所有图族均占优，也不能将原 24/32 线程的 2.902 倍移用到这里。

这些运行没有独立大规模 DP 核验；两方法共享 P3。
完整逐实例秒数、峰值、算法分支、来源类型和双方时间戳见 pairs.csv。
原始数据与来源表在 data/，源码在 code/，协议、环境与原验证记录在 provenance/。
'''
(HERE/'practical32'/'RESULTS.md').write_text(practical_md,encoding='utf-8')
verification=dict(input_files_verified=len(manifest),formal_pairs=len(theory_pairs),formal_runs=len(formal),
                  dp_runs=len(small),dp_graphs=len(small_pairs),practical_records=len(practical_records),
                  practical_pairs=len(practical_pairs),practical_graphs=len({(r['n'],r['graph_hash']) for r in practical_records}),
                  practical_shortcuts=sum(bool(r['shortcut']) for r in practical_pairs),
                  practical_sources=dict(Counter(r['source_mode'] for r in practical_pairs)),
                  all_archived_summary_values_match=True,all_checked_answers_and_hashes_match=True,
                  new_performance_runs=0)
write_json(HERE/'validation.json',verification)
print(json.dumps(verification,ensure_ascii=False))

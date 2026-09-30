# 两组独立对比材料

这些材料不属于论文正文。性能数字来自已有远程 i9-14900K 实验；本轮只复算数据与核验数学接口，没有新增性能测量。

## 理论配置

[结果说明](theory/RESULTS.md)：BH 与本项目均为单线程、完整确定化、多项式工作空间；压位 DP 是单线程确定性的指数空间参照。

| 入口 | 内容 |
|---|---|
| [summary.csv](theory/summary.csv) | 20–44 点规模汇总、输入范围、RSS、几何平均比值与胜数 |
| [pairs_dp.csv](theory/pairs_dp.csv) | 20–32 点三方法同轮重测，每图三次重复的中位数与范围 |
| [pairs_original.csv](theory/pairs_original.csv) | 原两方法 49 个正式配对及阶段时间 |
| [breadth.csv](theory/breadth.csv) | 36/40 点其他密度与 tournament |
| [时间与内存图](theory/figures/theory_comparison.pdf) | PDF；同目录另有 SVG、PNG |
| [data](theory/data) | 原始 JSONL 和原参考汇总 |
| [code](theory/code) | BH、本项目、压位 DP 及输入图生成源码 |
| [provenance](theory/provenance) | 协议、环境、原构建哈希与验证结果 |

20–32 点先按每图三次取中位数，再在三张图间取中位数；34–42 点每规模三张图各一次，44 点仅一图。
原两方法为 18 次训练加 98 次评估；新增三方对照为已有 21 张图上的 189 次，汇总没有重复计为新的独立输入。
时间比先配对后取几何平均，不等于两列中位数之比。图中误差棒是输入间最小到最大值，不是置信区间；DP 只画到实测 32 点。
计时实现的局部缓存保守为 `O(n³ log n)` 位，与论文严格流式 `O(n²)` 位实现区分。

[IMPLEMENTATION_NOTES.md](theory/IMPLEMENTATION_NOTES.md)保留确定化优化及空间依据；
[ADDITIONAL_THEORY.md](theory/ADDITIONAL_THEORY.md)保留原项目其他有用推导，作为参考材料，不是论文正文。

## 实际配置

[结果说明](practical32/RESULTS.md)：双方 32 硬件线程、Las Vegas、允许指数空间。
本项目保存访问列表以分桶精确去重；BH 实际仍用多项式空间，二部图使用其专用 Algorithm B。

| 入口 | 内容 |
|---|---|
| [summary.csv](practical32/summary.csv) | 规模、密度、结构、算法种子四组统计 |
| [pairs.csv](practical32/pairs.csv) | 138 个配对的秒数、内存、来源类型、时间戳和算法分支 |
| [时间与内存图](practical32/figures/practical32_comparison.pdf) | PDF；同目录另有 SVG、PNG |
| [data](practical32/data) | 原始组合数据、新增 BH32 记录及来源表 |
| [code](practical32/code) | 实际组源码、公共预处理、BH 二部专用实现 |
| [provenance](practical32/provenance) | 环境、原构建哈希、验证与完成记录 |
| [实现说明](practical32/IMPLEMENTATION_NOTES.md) | 随机空间与精确分桶去重的依据 |

138 个比较实例涉及 120 张不同图：133 组是补测 BH32 配已有本项目32记录，5 组复用已有完整 32/32 配对。
来源规则预先固定，不按时间择优。133 组双方测量时段不同，不能称为全组同时期交替重测；CPU 频率未锁定。
276 条组合记录不是 276 次新增实验。九个公共预处理捷径从核心速度比中排除，余下 129 组实际枚举。
大规模双方共享 P3，没有新增独立 DP 核验。50 点本项目峰值约 38.2 GiB，来自约 125 GiB 内存的远程机器。

## 复算、核验与文件来源

在 `v4` 内执行：

```powershell
python -B artifacts/build_artifacts.py
python -B artifacts/verify_math.py
```

第一个脚本读取两组 `data/`，核对 [inputs_manifest.json](inputs_manifest.json) 中的输入哈希，重算全部表图，
并逐项与原参考汇总核对。结果见 [validation.json](validation.json)。所有输出都在 v4 内。
第二个脚本只用 Python 标准库，核验全部 4164 张 2–4 点无自环有向图、18 张 5–7 点图、5421 个覆盖/归属组合、
24 个一般乘积系统的 720 次平移，以及控制器和精确矩恒等式。结果见 [mathematical_validation.json](mathematical_validation.json)。
这些是有限正确性检查，不是性能测量或渐近证明的替代物。

`inputs_manifest.json` 同时保留原文件路径、原 SHA-256、新路径与新 SHA-256。
理论组 C++ 只将 `../practical/baseline.cpp` 的 include 改为同目录 `baseline.cpp`，便于独立编译；算法内容未改。
原哈希记录中的旧文件名属于历史测量语境；本附件提供源码，没有附带原远程二进制。
原 `source` 目录及旧的依赖路径已移除，复算不再需要它。

如以后在足够内存的 Linux CPU 机器重现实验，可从 `artifacts` 目录编译到独立 `rebuilt` 目录：

```bash
mkdir -p rebuilt
g++ -O3 -march=native -std=c++17 -fopenmp theory/code/bench_poly.cpp -o rebuilt/bench_poly
g++ -O3 -march=native -std=c++17 -fopenmp theory/code/baseline.cpp -o rebuilt/bench_dp
g++ -O3 -march=native -std=c++17 -fopenmp practical32/code/extended.cpp -o rebuilt/extended
OMP_NUM_THREADS=1 OMP_DYNAMIC=FALSE OMP_PROC_BIND=TRUE OMP_PLACES='{0}' taskset -c 0 ./rebuilt/bench_poly kw-ce 26 9001 dense half
OMP_NUM_THREADS=1 OMP_DYNAMIC=FALSE OMP_PROC_BIND=TRUE OMP_PLACES='{0}' taskset -c 0 ./rebuilt/bench_poly bh-ce 26 9001 dense half
OMP_NUM_THREADS=1 OMP_DYNAMIC=FALSE OMP_PROC_BIND=TRUE OMP_PLACES='{0}' taskset -c 0 ./rebuilt/bench_dp dp-bit 26 9001 dense 1 two
```

实际组命令格式为 `./rebuilt/extended kw-practical 40 1001 dense 32 10001`，BH 对应 `bh-practical`。
完整种子与资源限制保存在数据和协议中；移机时按拓扑设置 CPU 绑定。这些命令本轮未运行。

"""Independent finite checks of the manuscript's mathematical interfaces.

Uses only Python's standard library. Results are not performance benchmarks.
The scalable C++ implementations are in theory/code/ and practical32/code/.
"""
from pathlib import Path
from itertools import product
from fractions import Fraction
from functools import lru_cache
import json
import random
import time

HERE=Path(__file__).resolve().parent
def parity(x): return x.bit_count()&1
def mask(n): return (1<<n)-1
def mv(rows,x): return sum(parity(row&x)<<i for i,row in enumerate(rows))
def solve(rows, rhs, n):
    a=[row|(((rhs>>i)&1)<<n) for i,row in enumerate(rows)]
    rank=0; piv=[]
    for col in range(n):
        at=next((j for j in range(rank,len(a)) if (a[j]>>col)&1),None)
        if at is None: continue
        a[rank],a[at]=a[at],a[rank]
        for j in range(len(a)):
            if j!=rank and ((a[j]>>col)&1): a[j]^=a[rank]
        piv.append(col);rank+=1
    if any((r&mask(n))==0 and ((r>>n)&1) for r in a): return None,rank
    return sum(((a[j]>>n)&1)<<c for j,c in enumerate(piv)),rank

def graph(n,bits,undirected=False):
    rows=[0]*n;k=0
    for i in range(n):
        for j in range(i+1 if undirected else 0,n):
            if i==j: continue
            if (bits>>k)&1:
                rows[i]|=1<<j
                if undirected:rows[j]|=1<<i
            k+=1
    return tuple(rows)

def centers(L):
    b=len(L);out=[]
    for x in range(1<<b):
        v=mv(L,x)
        if v&x==0:
            out.append(tuple(0 if x>>i&1 else 1 if v>>i&1 else 2 for i in range(b)))
    return out

def owner(L,s):
    rows=[];rhs=0
    for i in range(len(L)):
        rows.append(1<<i if s[i]==0 else L[i]^((1<<i) if s[i]==2 else 0))
        rhs|=(s[i]==2)<<i
    x,_=solve(rows,rhs,len(L));assert x is not None
    v=mv(L,x)
    assert x&v==0
    q=tuple(0 if x>>i&1 else 1 if v>>i&1 else 2 for i in range(len(L)))
    assert all(a!=b for a,b in zip(q,s))
    return q

@lru_cache(None)
def local_graph(b):
    # Exact completion averaging independently verifies the existence bound.
    candidates=[graph(b,x,True) for x in range(1<<(b*(b-1)//2))]
    L=min(candidates,key=lambda g:(len(centers(g)),g))
    assert len(centers(L))<=2*Fraction(3,2)**b-1
    return L

@lru_cache(None)
def cover(n):
    blocks=[];left=n
    while left:
        b=(left+1)//2;blocks.append(local_graph(b));left-=b
    locals_=[centers(g) for g in blocks]
    qs=[sum(parts,()) for parts in product(*locals_)]
    own={}
    for s in product(range(3),repeat=n):
        i=0;parts=[]
        for L in blocks:
            parts.append(owner(L,s[i:i+len(L)]));i+=len(L)
        own[s]=sum(parts,())
    return blocks,locals_,qs,own

class Rollback:
    def __init__(self,n):self.slots=[0]*n;self.stack=[]
    def residue(self,x):
        while x:
            p=x.bit_length()-1
            if not self.slots[p]:return x
            x^=self.slots[p]
        return 0
    def insert(self,x):
        x=self.residue(x)
        if x:
            p=x.bit_length()-1;self.slots[p]=x;self.stack.append(p)
    def undo(self,mark):
        while len(self.stack)>mark:self.slots[self.stack.pop()]=0

def prefix_counts(B,ell,c,blocks,lists):
    n=len(B);columns=[sum(((row>>i)&1)<<j for j,row in enumerate(B)) for i in range(n)]
    basis=Rollback(ell);counts=[0,0];insertions=0
    def visit(j,start,shift):
        nonlocal insertions
        if j==len(blocks):
            for bit in [0,1]:
                rhs=(c|(bit<<(ell-1)))^shift
                if basis.residue(rhs)==0:counts[bit]+=1<<(n-len(basis.stack))
            return
        for q in lists[j]:
            mark=len(basis.stack);d=shift
            for k,sym in enumerate(q):
                i=start+k;v=columns[i]^(1<<i) if sym==0 else (1<<i) if sym==1 else columns[i]
                basis.insert(v&mask(ell));insertions+=1
                if sym==0:d^=columns[i]&mask(ell)
            visit(j+1,start+len(q),d)
            basis.undo(mark)
    visit(0,0,0)
    return counts,insertions

def p2(B,c):return [x for x in range(1<<len(B)) if x&(mv(B,x)^c)==0]

def weight(B,c,x):
    n=len(B);A=[row^((((c>>i)&1)^1)<<i) for i,row in enumerate(B)]
    low=(x&-x).bit_length()-1
    rows=[];rhs=[]
    for i in range(n):
        if (x>>i)&1 or i<=low:rows.append(1<<i);rhs.append(0)
        if not (x>>i)&1:
            d=parity(A[i]&x);rows.append(A[i]^(d<<i));rhs.append(1^d)
    sol,rank=solve(rows,sum(v<<i for i,v in enumerate(rhs)),n)
    return int(sol is not None and rank==n)

def subset_dp(B):
    n=len(B);D={1:1}
    for s in range(1,1<<n,2):
        bits=D.get(s,0)
        for v in range(1,n):
            if s>>v&1:continue
            bit=0
            for u in range(n):bit^=((bits>>u)&1)&((B[u]>>v)&1)
            if bit:D[s|(1<<v)]=D.get(s|(1<<v),0)^(1<<v)
    return parity(D.get(mask(n),0)&sum(((B[v]&1)>0)<<v for v in range(n)))

def verify_graph(B):
    n=len(B);blocks,lists,qs,owners=cover(n)
    columns=[sum(((row>>i)&1)<<j for j,row in enumerate(B)) for i in range(n)]
    # Brute-force all (center, parameter) pairs, independently of elimination.
    profile=[0]*(1<<n);visits=[[] for _ in range(1<<n)]
    for q in qs:
        for z in range(1<<n):
            x=y=0
            for i,sym in enumerate(q):
                bit=(z>>i)&1
                if sym==0:x|=(1^bit)<<i;y|=bit<<i
                elif sym==1:y|=bit<<i
                else:x|=bit<<i
            syndrome=mv(B,x)^y
            profile[syndrome]+=1;visits[syndrome].append((q,x,y))
    assert sum(profile)==len(qs)*(1<<n)
    c=0;maxratio=0
    for ell in range(1,n+1):
        expected=[sum(profile[d] for d in range(1<<n) if (d&mask(ell))==(c|(b<<(ell-1)))) for b in [0,1]]
        actual,I=prefix_counts(B,ell,c,blocks,lists)
        assert actual==expected,(B,ell,actual,expected)
        maxratio=max(maxratio,Fraction(I,len(qs)));assert I<=9*len(qs)
        if actual[1]<actual[0]:c|=1<<(ell-1)
    assert profile[c]<=len(qs)
    kept=[]
    for q,x,y in visits[c]:
        s=tuple(1 if x>>i&1 else 2 if y>>i&1 else 0 for i in range(n))
        if q==owners[s]:kept.append(x)
    assert sorted(kept)==p2(B,c)
    result=0
    for x in kept:
        if x:result^=weight(B,c,x)
    assert result==subset_dp(B),(B,c,result)
    # Every diagonal produces the same weighted answer on the smallest graphs.
    if n<=3:
        for d in range(1<<n):
            ans=0
            for x in p2(B,d):
                if x:ans^=weight(B,d,x)
            assert ans==result
    return maxratio

def main():
    started=time.perf_counter();ownership_checks=0
    for b in range(1,5):
        for bits in range(1<<(b*(b-1)//2)):
            L=graph(b,bits,True);qs=centers(L)
            for s in product(range(3),repeat=b):
                assert owner(L,s) in qs;ownership_checks+=1
    exhaustive=0;maxratio=Fraction(0)
    for n in range(2,5):
        for bits in range(1<<(n*(n-1))):
            maxratio=max(maxratio,verify_graph(graph(n,bits)));exhaustive+=1
        print(f'Checked all loopless digraphs at n={n}.',flush=True)
    rng=random.Random(159204);sampled=0
    for n in [5,6,7]:
        for _ in range(6):
            B=tuple(rng.getrandbits(n)&~(1<<i) for i in range(n))
            maxratio=max(maxratio,verify_graph(B));sampled+=1
    # Exact joint first and second moments on all n<=3 matrices and syndromes.
    for n in [2,3]:
        counts=[len(p2(graph(n,bits),c)) for bits in range(1<<(n*(n-1))) for c in range(1<<n)]
        mean=Fraction(sum(counts),len(counts));var=sum((v-mean)**2 for v in counts)/len(counts)
        assert mean==Fraction(3,2)**n
        assert var==Fraction(3,2)**n-Fraction(5,4)**n
    G=(6,5,3,48,40,24)
    controller_counts=[len(p2(G,c)) for c in range(64)]
    assert min(controller_counts)==9
    assert 48*8**32>=9**32 and 48*8**33<9**33
    assert 1008*8**58>=9**58 and 1008*8**59<9**59
    assert 8*Fraction(2,3)**7==Fraction(1024,2187)
    systems=24;translations=0
    for case in range(systems):
        N=3+case%3;m=1+case%4
        forms=[(rng.getrandbits(N),rng.getrandbits(1),rng.getrandbits(N),rng.getrandbits(1)) for _ in range(m)]
        def state(x):
            pairs=[(parity(u&x)^uc,parity(v&x)^vc) for u,uc,v,vc in forms]
            if any(a&b for a,b in pairs):return None
            return tuple(1 if a else 2 if b else 0 for a,b in pairs)
        states={x:state(x) for x in range(1<<N)};legal={x:s for x,s in states.items() if s is not None}
        _,_,qs,owners=cover(m);total_visits=0
        for t in product(range(3),repeat=m):
            outputs=[]
            for q in qs:
                qt=tuple((a+b)%3 for a,b in zip(q,t))
                for x,s in legal.items():
                    if all(a!=b for a,b in zip(qt,s)):
                        total_visits+=1
                        if owners[tuple((a-b)%3 for a,b in zip(s,t))]==q:outputs.append(x)
            assert sorted(outputs)==sorted(legal)
            translations+=1
        assert total_visits==len(qs)*(2**m)*len(legal)
    result=dict(status='passed',exhaustive_loopless_graphs=exhaustive,additional_graphs=sampled,
                ownership_graph_state_pairs=ownership_checks,max_insertion_to_center_ratio=str(maxratio),
                product_systems=systems,translations=translations,controller_syndromes=64,
                minimum_controller_solutions=min(controller_counts),exact_thresholds=[33,59],
                checks=['rollback prefix counts versus exhaustive parameter enumeration',
                        'visit budget, ownership uniqueness, weighted P3 versus independent subset DP',
                        'all diagonals at n<=3', 'exact joint moments', 'controller amplification constants',
                        'general product listing and translation visit identity'],
                seconds=round(time.perf_counter()-started,3),performance_benchmark=False)
    (HERE/'mathematical_validation.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(result),flush=True)

if __name__=='__main__':main()

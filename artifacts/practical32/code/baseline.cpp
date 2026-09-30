// CPU-only research benchmark, C++17/OpenMP, n <= 40.
// Formulas: BH13 Algorithms C/D and P3; KW26 Lemma 2.2; 1592/v3 and v3-ex1..4.
// All methods compute the complete Hamiltonian-cycle parity, not P2 parity.
#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>
#include <sys/resource.h>
#include <omp.h>
using U=uint64_t;
using W=unsigned __int128;
using Clock=std::chrono::steady_clock;
static int pc(U x){return __builtin_popcountll(x);}
static int par(U x){return __builtin_parityll(x);}
static int low(U x){return __builtin_ctzll(x);}
static U mask(int n){return (U(1)<<n)-1;}
static double secs(Clock::time_point t){return std::chrono::duration<double>(Clock::now()-t).count();}
struct Graph{int n; std::vector<U> b,col;};
static Graph graph(int n,U seed,std::string kind){
    Graph g{n,std::vector<U>(n),std::vector<U>(n)};
    std::mt19937_64 rng(seed);
    for(int i=0;i<n;i++)for(int j=0;j<n;j++)if(i!=j){
        bool edge=false;
        U r=rng();
        if(kind=="dense")edge=(r%1000)<500;
        else if(kind=="sparse")edge=(r%1000)<150;
        else if(kind=="verydense")edge=(r%1000)<850;
        else if(kind=="cycle")edge=j==(i+1)%n;
        else if(kind=="complete")edge=true;
        else if(kind=="bipartite")edge=((i<n/2)!=(j<n/2))&&(r%1000<500);
        else if(kind=="dag")edge=j>i && (r%1000<500);
        else if(kind=="tournament"){if(i<j){if(r&1)g.b[i]|=U(1)<<j;else g.b[j]|=U(1)<<i;}continue;}
        else if(kind=="controller"){
            int m=n-7;
            if(i<m)edge=(r&1);
            else if(i<n-1 && j>=m && j<n-1)edge=(i-m)/3==(j-m)/3;
            if(i==n-1 || j==n-1)edge=true;
        }else if(kind!="empty")throw std::runtime_error("unknown graph kind");
        if(edge)g.b[i]|=U(1)<<j;
    }
    for(int i=0;i<n;i++)for(int j=0;j<n;j++)if(g.b[i]>>j&1)g.col[j]|=U(1)<<i;
    return g;
}
static void columns(Graph& g){g.col.assign(g.n,0);for(int i=0;i<g.n;i++)for(int j=0;j<g.n;j++)if(g.b[i]>>j&1)g.col[j]|=U(1)<<i;}
struct Basis{
    U row[64]{}; int rank=0;
    // Augmented row: RHS in bit n. Only insertions alter the basis.
    int add(U a,int rhs,int n){
        U v=a|(U(rhs)<<n);
        while(v&mask(n)){int p=low(v);if(!row[p]){row[p]=v;rank++;return p;}v^=row[p];}
        return v?-2:-1;
    }
    U reduce(U a,int rhs,int n)const{
        U v=a|(U(rhs)<<n);
        for(int p=0;p<n;p++)if((v>>p&1)&&row[p])v^=row[p];
        return v;
    }
    void undo(int p){if(p>=0){row[p]=0;rank--;}}
    U particular(int n)const{
        U x=0;
        for(int p=n-1;p>=0;p--)if(row[p])if(par(row[p]&x)^int(row[p]>>n&1))x|=U(1)<<p;
        return x;
    }
    template<class F>void stream(int n,F f)const{
        U x=particular(n), dirs[64];int d=0;
        for(int j=0;j<n;j++)if(!row[j]){
            U z=U(1)<<j;
            for(int p=n-1;p>=0;p--)if(row[p] && par(row[p]&z))z|=U(1)<<p;
            dirs[d++]=z;
        }
        U total=U(1)<<d;
        for(U i=0;i<total;i++){if(i)x^=dirs[low(i)];f(x);}
    }
};
static bool p3(const Graph&g,U c,U x){
    if(!x)return false;
    int n=g.n, first=low(x);U free=mask(n)&~(x|mask(first+1));
    Basis b;
    U rest=mask(n)^x;
    while(rest){int i=low(rest);rest&=rest-1;int d=par(g.b[i]&x);
        U a=(g.b[i]^(((c>>i&1)^1^d)<<i))&free;
        if(b.add(a,1^d,n)==-2)return false;
    }
    return b.rank==pc(free);
}
static bool legal(const Graph&g,U c,U x){U r=x;while(r){int i=low(r);r&=r-1;if(par(g.b[i]&x)!=int(c>>i&1))return false;}return true;}
struct Result{
    int answer=0;U c=0,centers=0,visits=0,unique=0,nodes=0,tablebytes=0;
    double setup=0,choose=0,enumerate=0,dedup=0,weight=0,total=0;
};
static Result dp_naive(const Graph&g){
    Result r;int n=g.n,m=n-1;U count=U(1)<<m;
    std::vector<unsigned char>d(size_t(count)*m);
    for(int v=0;v<m;v++)d[(U(1)<<v)*m+v]=g.b[0]>>(v+1)&1;
    for(U s=1;s<count;s++)for(int v=0;v<m;v++)if((s>>v&1)&&d[s*m+v]){
        for(int w=0;w<m;w++)if(!(s>>w&1)&&(g.b[v+1]>>(w+1)&1))d[(s|(U(1)<<w))*m+w]^=1;
    }
    for(int v=0;v<m;v++)if(g.b[v+1]&1)r.answer^=d[(count-1)*m+v];
    r.tablebytes=count*m;return r;
}
static Result dp_bit(const Graph&g,int threads){
    Result r;int m=g.n-1;U count=U(1)<<m;
    std::vector<U>d(count);U in[64]{};
    for(int v=0;v<m;v++)in[v]=g.col[v+1]>>1;
    auto process=[&](U s){U bits=s,val=0;while(bits){int v=low(bits);bits&=bits-1;U prev=s^(U(1)<<v);int x=prev?par(d[prev]&in[v]):int(g.b[0]>>(v+1)&1);val|=U(x)<<v;}d[s]=val;};
    if(threads==1){for(U s=1;s<count;s++)process(s);}
    else{for(int k=1;k<=m;k++){
        #pragma omp parallel for schedule(static)
        for(U s=1;s<count;s++)if(pc(s)==k)process(s);
    }}
    r.answer=par(d.back()&(g.col[0]>>1));r.tablebytes=count*8;return r;
}
static U binom[64][64];
static void init_binom(){for(int n=0;n<64;n++){binom[n][0]=binom[n][n]=1;for(int k=1;k<n;k++)binom[n][k]=binom[n-1][k-1]+binom[n-1][k];}}
static U unrank_comb(int n,int k,U r){U s=0;for(int j=k;j>=1;j--){int p=n-1;while(binom[p][j]>r)p--;s|=U(1)<<p;r-=binom[p][j];n=p;}return s;}
static U nextcomb(U x){U u=x&-x,v=x+u;return v+(((v^x)/u)>>2);}
struct Prefix{U x;int len,k;U id=0;};
template<class F>static void bh_jobs(int n,F f){
    U offset=0;
    for(int k=0;k<=n/2;k++)for(int type=0;type<2;type++){
        int len=n-k,positions=len-type;if(positions<k)continue;U count=binom[positions][k];
        #pragma omp parallel for schedule(dynamic,1)
        for(U base=0;base<count;base+=256){U s=unrank_comb(positions,k,base);for(U z=base;z<std::min(base+256,count);z++){
            f(Prefix{s|(type?(U(1)<<(len-1)):0),len,k,offset+z});if(k && z+1<std::min(base+256,count))s=nextcomb(s);
        }}
        offset+=count;
    }
}
static W bh_count(const Graph&g,Prefix p,U c,int ell){
    Basis b;U active=p.x;int fixed=0;
    while(active){int i=low(active);active&=active-1;if(i>=ell)continue;fixed++;
        U a=g.b[i]>>p.len;int rhs=par(g.b[i]&p.x)^int(c>>i&1);
        if(b.add(a,rhs,p.k)==-2)return 0;
    }
    int unknown=pc(p.x)-fixed;
    return W(1)<<(g.n-ell+p.k-b.rank-unknown);
}
static Result bh(const Graph&g,U seed,bool deterministic,bool forced=false,U forced_c=0){
    Result r;int n=g.n;std::mt19937_64 rng(seed^0xB113);r.c=rng()&mask(n);
    if(forced)r.c=forced_c;
    auto t=Clock::now();
    if(deterministic){r.c=0;for(int ell=1;ell<=n;ell++){
        struct alignas(64) Sum{W a=0,b=0;};std::vector<Sum>sums(omp_get_max_threads());
        bh_jobs(n,[&](Prefix p){auto&s=sums[omp_get_thread_num()];s.a+=bh_count(g,p,r.c,ell);s.b+=bh_count(g,p,r.c|(U(1)<<(ell-1)),ell);});
        W a=0,b=0;for(auto s:sums){a+=s.a;b+=s.b;}if(b<a)r.c|=U(1)<<(ell-1);
    }}r.choose=secs(t);t=Clock::now();
    struct alignas(64) Local{U centers=0,visits=0,unique=0;int h=0;};
    std::vector<Local> locals(omp_get_max_threads());
    bh_jobs(n,[&](Prefix p){auto&l=locals[omp_get_thread_num()];l.centers++;Basis b;U active=p.x;
        while(active){int i=low(active);active&=active-1;if(b.add(g.b[i]>>p.len,par(g.b[i]&p.x)^int(r.c>>i&1),p.k)==-2)return;}
        b.stream(p.k,[&](U tail){l.visits++;U x=p.x|(tail<<p.len);if(legal(g,r.c,x)){l.unique++;l.h^=p3(g,r.c,x);}});
    });
    for(auto l:locals){r.centers+=l.centers;r.visits+=l.visits;r.unique+=l.unique;r.answer^=l.h;}
    r.enumerate=secs(t);return r;
}
static Result bh_cached(const Graph&g,U seed){
    auto t=Clock::now();int n=g.n,threads=omp_get_max_threads();U M=0;
    for(int k=0;k<=n/2;k++)for(int type=0;type<2;type++)if(n-k-type>=k)M+=binom[n-k-type][k];
    struct Rec{U id_depth,affine;};
    std::vector<std::vector<std::vector<Rec>>>local(threads,std::vector<std::vector<Rec>>(n));
    bh_jobs(n,[&](Prefix p){
        U rows[64]{},tags[64]{};int right[64]{};int depth=1+p.k-pc(p.x);U bits=p.x;
        while(bits){int i=low(bits);bits&=bits-1;U a=g.b[i]>>p.len,tag=U(1)<<i;int rhs=par(g.b[i]&p.x);
            while(a){int pivot=low(a);if(!rows[pivot]){rows[pivot]=a;tags[pivot]=tag;right[pivot]=rhs;break;}a^=rows[pivot];tag^=tags[pivot];rhs^=right[pivot];}
            if(!a)local[omp_get_thread_num()][i].push_back({(p.id<<6)|U(depth++),tag|(U(rhs)<<n)});
        }
    });
    U c=0,count=0;
    {
        std::vector<std::vector<Rec>>deps(n);
        for(int i=0;i<n;i++){U size=0;for(auto&t:local)size+=t[i].size();deps[i].reserve(size);for(auto&t:local){deps[i].insert(deps[i].end(),t[i].begin(),t[i].end());std::vector<Rec>().swap(t[i]);}count+=size;}
        std::vector<uint8_t>alive(M,1);
        for(int i=0;i<n;i++){
            __int128 delta=0;
            #pragma omp parallel for reduction(+:delta) schedule(static)
            for(U j=0;j<deps[i].size();j++){auto v=deps[i][j];if(alive[v.id_depth>>6]){int req=par(c&v.affine)^int(v.affine>>n&1);__int128 w=__int128(1)<<(v.id_depth&63);delta+=req?-w:w;}}
            if(delta>0)c|=U(1)<<i;
            #pragma omp parallel for schedule(static)
            for(U j=0;j<deps[i].size();j++){auto v=deps[i][j];if(alive[v.id_depth>>6]&&(par(c&v.affine)^int(v.affine>>n&1)))alive[v.id_depth>>6]=0;}
        }
    }
    double choose=secs(t);auto r=bh(g,seed,false,true,c);r.choose=choose;r.tablebytes=M+48*count;return r;
}
struct Center{U z=0,o=0;}; // prohibited 0 and prohibited 1; other coordinates prohibit 2
struct Block{int b,off;std::vector<U> graph;std::vector<Center> list;U count=0;std::vector<int> owners;};
static Center center_for(const Block&b,U s){Center q{s,0};for(int i=0;i<b.b;i++)if(!(s>>i&1)&&par(b.graph[i]&s))q.o|=U(1)<<i;return q;}
static bool even(const Block&b,U s){U bits=s;while(bits){int i=low(bits);bits&=bits-1;if(par(b.graph[i]&s))return false;}return true;}
// Exact edge conditional expectation: scan subsets, connected components of unexposed edges.
static W graph_ce(int n,const std::vector<U>&fixed,const std::vector<U>&unknown){
    W total=0;
    for(U s=0;s<(U(1)<<n);s++){
        U rem=s;int comps=0;bool ok=true;
        while(rem){U component=U(1)<<low(rem),front=component;rem^=component;
            while(front){int i=low(front);front&=front-1;U add=unknown[i]&rem;component|=add;front|=add;rem^=add;}
            int h=0;U bits=component;while(bits){int i=low(bits);bits&=bits-1;h^=par(fixed[i]&s);}if(h){ok=false;break;}comps++;
        }if(ok)total+=W(1)<<(n-pc(s)+comps);
    }return total;
}
static Block build_block(int b,int off,bool det,U seed,bool cache=true){
    Block out{b,off,std::vector<U>(b),{},0,{}};
    if(det){std::vector<U> unknown(b);for(int i=0;i<b;i++)unknown[i]=mask(b)^(U(1)<<i);
        for(int i=0;i<b;i++)for(int j=i+1;j<b;j++){
            unknown[i]^=U(1)<<j;unknown[j]^=U(1)<<i;
            W s0=graph_ce(b,out.graph,unknown);
            out.graph[i]^=U(1)<<j;out.graph[j]^=U(1)<<i;
            W s1=graph_ce(b,out.graph,unknown);
            if(s0<=s1){out.graph[i]^=U(1)<<j;out.graph[j]^=U(1)<<i;}
        }
    }else{
        // One graph, no uncharged best-of-many auxiliary search.
        std::mt19937_64 rng(seed);for(int i=0;i<b;i++)for(int j=i+1;j<b;j++)if(rng()&1){out.graph[i]|=U(1)<<j;out.graph[j]|=U(1)<<i;}
    }
    for(U s=0;s<(U(1)<<b);s++)if(even(out,s)){out.count++;if(cache)out.list.push_back(center_for(out,s));}
    return out;
}
static Center global(const std::vector<Block>&bs,U index){Center q;for(int j=int(bs.size())-1;j>=0;j--){auto&b=bs[j];auto c=b.list[index%b.list.size()];index/=b.list.size();q.z|=c.z<<b.off;q.o|=c.o<<b.off;}return q;}
// Generic greedy cover: scores=(J_3-I_3)^{tensor b} uncovered.
// O(b^2 (9/2)^b) time, O(b 3^b) bits. Four balanced blocks keep setup below 1.5^n.
static Block greedy_block(int b,int off){
    Block out{b,off,{},{},0,{}};int power[32],N=1;for(int i=0;i<b;i++){power[i]=N;N*=3;}
    out.owners.assign(N,-1);std::vector<int>uncovered(N,1),scores(N);int remaining=N;
    while(remaining){
        scores=uncovered;
        for(int stride=1;stride<N;stride*=3)for(int base=0;base<N;base+=3*stride)for(int j=0;j<stride;j++){
            int a=scores[base+j],c=scores[base+j+stride],d=scores[base+j+2*stride];
            scores[base+j]=c+d;scores[base+j+stride]=a+d;scores[base+j+2*stride]=a+c;
        }
        int best=int(std::max_element(scores.begin(),scores.end())-scores.begin()),v=best,index=0,delta[32];Center q;
        for(int i=0;i<b;i++){int s=v%3;v/=3;if(s==0)q.z|=U(1)<<i;if(s==1)q.o|=U(1)<<i;
            int a=s==0?1:0,c=s==2?1:2;index+=a*power[i];delta[i]=(c-a)*power[i];
        }
        int owner=int(out.list.size());out.list.push_back(q);U gray=0;
        for(U step=0;step<(U(1)<<b);step++){
            if(step){int bit=low(step);gray^=U(1)<<bit;index+=(gray>>bit&1)?delta[bit]:-delta[bit];}
            if(uncovered[index]){uncovered[index]=0;out.owners[index]=owner;remaining--;}
        }
    }
    out.count=out.list.size();return out;
}
static U syndrome(const Graph&g,U x){U c=0;while(x){int i=low(x);x&=x-1;c^=g.col[i];}return c;}
static U qrow(const Graph&g,Center q,int i){return (g.b[i]&~q.o)^(((q.z|q.o)>>i&1)<<i);}
static bool qsolve(const Graph&g,Center q,U c,Basis&b){U rhs=c^syndrome(g,q.z);for(int i=0;i<g.n;i++)if(b.add(qrow(g,q,i),int(rhs>>i&1),g.n)==-2)return false;return true;}
static bool owned(const std::vector<Block>&bs,Center q,U x,U y){
    for(auto&b:bs){U xx=(x>>b.off)&mask(b.b), yy=(y>>b.off)&mask(b.b);Basis owner;
        if(!b.owners.empty()){
            int index=0,power=1;for(int i=0;i<b.b;i++){index+=power*((xx>>i&1)?1:((yy>>i&1)?2:0));power*=3;}
            auto o=b.list[b.owners[index]];if(o.z!=((q.z>>b.off)&mask(b.b))||o.o!=((q.o>>b.off)&mask(b.b)))return false;
            continue;
        }
        for(int i=0;i<b.b;i++){
            U a;int rhs;
            if(xx>>i&1){a=b.graph[i];rhs=0;}
            else if(yy>>i&1){a=b.graph[i]^(U(1)<<i);rhs=1;}
            else{a=U(1)<<i;rhs=0;}
            if(owner.add(a,rhs,b.b)==-2)throw std::runtime_error("owner inconsistent");
        }
        U s=owner.particular(b.b);if(s!=((q.z>>b.off)&mask(b.b)))return false;
    }return true;
}
static U choose_naive(const Graph&g,const std::vector<Block>&bs,U M){
    U c=0;for(int ell=1;ell<=g.n;ell++){W a=0,b=0;
        #pragma omp parallel for reduction(+:a,b) schedule(static)
        for(U id=0;id<M;id++){Center q=global(bs,id);U d=syndrome(g,q.z);
            for(int bit=0;bit<2;bit++){Basis base;U rhs=c^(U(bit)<<(ell-1))^d;bool ok=true;
                for(int i=0;i<ell;i++)if(base.add(qrow(g,q,i),int(rhs>>i&1),g.n)==-2){ok=false;break;}
                if(ok){W w=W(1)<<(g.n-base.rank);if(bit)b+=w;else a+=w;}
            }
        }if(b<a)c|=U(1)<<(ell-1);
    }return c;
}
static U choose_cache(const Graph&g,const std::vector<Block>&bs,U M,U&bytes){
    int n=g.n;std::vector<U> deps(size_t(M)*n),d(M);std::vector<uint8_t> rank(M),alive(M,1);
    bytes=U(deps.size()+d.size())*8+2*M;
    #pragma omp parallel for schedule(static)
    for(U id=0;id<M;id++){
        Center q=global(bs,id);d[id]=syndrome(g,q.z);U rows[64]{},tags[64]{};
        for(int i=0;i<n;i++){U a=qrow(g,q,i),tag=U(1)<<i;while(a){int p=low(a);if(!rows[p]){rows[p]=a;tags[p]=tag;break;}a^=rows[p];tag^=tags[p];}if(!a)deps[U(i)*M+id]=tag;}
    }
    U c=0;
    for(int i=0;i<n;i++){W a=0,b=0;
        #pragma omp parallel for reduction(+:a,b) schedule(static)
        for(U id=0;id<M;id++)if(alive[id]){
            U dep=deps[U(i)*M+id];if(!dep){rank[id]++;W w=W(1)<<(n-rank[id]);a+=w;b+=w;}
            else{W w=W(1)<<(n-rank[id]);if(par((c^d[id])&dep))b+=w;else a+=w;}
        }
        if(b<a)c|=U(1)<<i;
        #pragma omp parallel for schedule(static)
        for(U id=0;id<M;id++)if(alive[id]&&deps[U(i)*M+id]&&par((c^d[id])&deps[U(i)*M+id]))alive[id]=0;
    }return c;
}
// Normalize each prefix contribution to 2^(ell-rank), so an independent row
// changes nothing. Only dependent rows affect the two conditional expectations.
// A record saves center id, number of earlier dependencies, and one affine test.
static U choose_sparse(const Graph&g,const std::vector<Block>&bs,U M,U&bytes){
    struct Rec{U id_depth,affine;};int n=g.n,threads=omp_get_max_threads();
    std::vector<std::vector<std::vector<Rec>>>local(threads,std::vector<std::vector<Rec>>(n));
    #pragma omp parallel for schedule(static)
    for(U id=0;id<M;id++){
        Center q=global(bs,id);U d=syndrome(g,q.z),rows[64]{},tags[64]{};int depth=0;
        for(int i=0;i<n;i++){U a=qrow(g,q,i),tag=U(1)<<i;while(a){int p=low(a);if(!rows[p]){rows[p]=a;tags[p]=tag;break;}a^=rows[p];tag^=tags[p];}
            if(!a)local[omp_get_thread_num()][i].push_back({(id<<6)|U(depth++),tag|(U(par(d&tag))<<n)});
        }
    }
    std::vector<std::vector<Rec>>deps(n);U count=0;
    for(int i=0;i<n;i++){U size=0;for(auto&t:local)size+=t[i].size();deps[i].reserve(size);for(auto&t:local){deps[i].insert(deps[i].end(),t[i].begin(),t[i].end());std::vector<Rec>().swap(t[i]);}count+=size;}
    bytes=M+48*count; // Includes vector capacity slack and simultaneous merging storage.
    std::vector<uint8_t>alive(M,1);U c=0;
    for(int i=0;i<n;i++){
        __int128 delta=0;
        #pragma omp parallel for reduction(+:delta) schedule(static)
        for(U j=0;j<deps[i].size();j++){
            Rec v=deps[i][j];U id=v.id_depth>>6;if(!alive[id])continue;
            int required=par(c&v.affine)^int(v.affine>>n&1);__int128 weight=__int128(1)<<(v.id_depth&63);
            delta+=required?-weight:weight;
        }
        if(delta>0)c|=U(1)<<i;
        #pragma omp parallel for schedule(static)
        for(U j=0;j<deps[i].size();j++){
            Rec v=deps[i][j];U id=v.id_depth>>6;if(alive[id]&&(par(c&v.affine)^int(v.affine>>n&1)))alive[id]=0;
        }
    }return c;
}
static U choose_rollback(const Graph&g,const std::vector<Block>&bs,bool stream){
    int n=g.n;U c=0;
    for(int ell=1;ell<=n;ell++){
        W a=0,b=0;U bas[64]{};int rank=0;
        auto contains=[&](U x){while(x){int p=low(x);if(!bas[p])return false;x^=bas[p];}return true;};
        std::function<void(int,U)>dfs=[&](int j,U d){
            if(j==int(bs.size())){W w=W(1)<<(n-rank);if(contains(c^d))a+=w;if(contains(c^d^(U(1)<<(ell-1))))b+=w;return;}
            const auto&block=bs[j];
            auto visit=[&](Center q){int piv[64],cnt=0;U dd=d;
                for(int k=0;k<block.b;k++){int i=k+block.off;U v;
                    if(q.z>>k&1){v=g.col[i]^(U(1)<<i);dd^=g.col[i]&mask(ell);}
                    else if(q.o>>k&1)v=U(1)<<i;else v=g.col[i];
                    v&=mask(ell);
                    while(v){int p=low(v);if(!bas[p]){bas[p]=v;rank++;piv[cnt++]=p;break;}v^=bas[p];}
                }dfs(j+1,dd);while(cnt){bas[piv[--cnt]]=0;rank--;}
            };
            if(stream){for(U s=0;s<(U(1)<<block.b);s++)if(even(block,s))visit(center_for(block,s));}
            else for(Center q:block.list)visit(q);
        };dfs(0,0);if(b<a)c|=U(1)<<(ell-1);
    }return c;
}
static Result kw(const Graph&g,U seed,std::string method,std::string layout){
    Result r;int n=g.n;bool det=method.find("det")!=std::string::npos || method=="kw-naive" || method=="kw-roll" || method=="kw-stream";
    bool bucket=method.find("bucket")!=std::string::npos;
    bool sort=method.find("sort")!=std::string::npos||bucket,generic=method.rfind("greedy-",0)==0;
    auto t=Clock::now();std::vector<int> sizes;
    if(generic){for(int j=0;j<4;j++)sizes.push_back(n/4+(j<n%4));}
    else if(layout=="two")sizes={(n+1)/2,n/2};
    else if(layout=="geo"){int left=n;while(left){sizes.push_back((left+1)/2);left/=2;}}
    else if(layout=="six"){int left=n;while(left){int k=std::min(6,left);sizes.push_back(k);left-=k;}}
    else throw std::runtime_error("layout");
    std::vector<Block>bs;int off=0;r.centers=1;
    for(int sz:sizes)if(sz){
        int reuse=-1;
        if(method.find("sparse-det")!=std::string::npos)for(int j=0;j<int(bs.size());j++)if(bs[j].b==sz){reuse=j;break;}
        if(reuse>=0){Block copy=bs[reuse];copy.off=off;bs.push_back(std::move(copy));}
        else bs.push_back(generic?greedy_block(sz,off):build_block(sz,off,det,seed^U(sz*1009+off),method!="kw-stream"));
        off+=sz;r.centers*=bs.back().count;
    }
    r.setup=secs(t);t=Clock::now();
    if(method=="kw-naive"||method=="greedy-det-lowmem")r.c=choose_naive(g,bs,r.centers);
    else if(method=="kw-roll"||method=="kw-stream")r.c=choose_rollback(g,bs,method=="kw-stream");
    else if(det && method.find("sparse")!=std::string::npos)r.c=choose_sparse(g,bs,r.centers,r.tablebytes);
    else if(det)r.c=choose_cache(g,bs,r.centers,r.tablebytes);
    else{std::mt19937_64 rng(seed^0xB113);r.c=rng()&mask(n);}
    r.choose=secs(t);t=Clock::now();
    struct alignas(64) Local{U visits=0,unique=0;int h=0;std::vector<U> xs;};std::vector<Local>locals(omp_get_max_threads());
    auto process=[&](Center q){
        auto&l=locals[omp_get_thread_num()];Basis b;if(!qsolve(g,q,r.c,b))return;
        b.stream(n,[&](U z){U x=(z^q.z)&~q.o;l.visits++;
            if(sort)l.xs.push_back(x);
            else{U y=z&(q.z|q.o);if(owned(bs,q,x,y)){l.unique++;l.h^=p3(g,r.c,x);}}
        });
    };
    if(method=="kw-stream"){
        std::function<void(int,Center)>dfs=[&](int j,Center q){
            if(j==int(bs.size())){process(q);return;}
            auto&b=bs[j];for(U s=0;s<(U(1)<<b.b);s++)if(even(b,s)){auto local=center_for(b,s);dfs(j+1,{q.z|(local.z<<b.off),q.o|(local.o<<b.off)});}
        };dfs(0,{});
    }else{
        #pragma omp parallel for schedule(dynamic,128)
        for(U id=0;id<r.centers;id++)process(global(bs,id));
    }
    r.enumerate=secs(t);for(auto&l:locals){r.visits+=l.visits;r.unique+=l.unique;r.answer^=l.h;}
    if(sort){t=Clock::now();std::vector<U>xs;xs.reserve(r.visits);for(auto&l:locals)xs.insert(xs.end(),l.xs.begin(),l.xs.end());
        U listbytes=24*r.visits;
        if(bucket)listbytes=std::max(listbytes,16*r.visits+U(omp_get_max_threads())*1024*16);
        r.tablebytes=std::max(r.tablebytes,listbytes);for(auto&l:locals){std::vector<U>().swap(l.xs);}
        if(bucket){
            constexpr int K=1024;int threads=omp_get_max_threads();
            std::vector<std::array<U,K>>hist(threads),offsets(threads);std::array<U,K+1>starts{};
            auto bin=[](U x){x^=x>>30;x*=0xbf58476d1ce4e5b9ULL;x^=x>>27;x*=0x94d049bb133111ebULL;x^=x>>31;return int(x&(K-1));};
            std::vector<U>scratch(xs.size());
            #pragma omp parallel
            {int tid=omp_get_thread_num();
                #pragma omp for schedule(static)
                for(U i=0;i<xs.size();i++)hist[tid][bin(xs[i])]++;
                #pragma omp single
                for(int k=0;k<K;k++){U count=0;for(int worker=0;worker<threads;worker++){offsets[worker][k]=starts[k]+count;count+=hist[worker][k];}starts[k+1]=starts[k]+count;}
                #pragma omp for schedule(static)
                for(U i=0;i<xs.size();i++)scratch[offsets[tid][bin(xs[i])]++]=xs[i];
            }
            xs.swap(scratch);std::vector<U>().swap(scratch);
            #pragma omp parallel for schedule(dynamic,1)
            for(int k=0;k<K;k++)std::sort(xs.begin()+starts[k],xs.begin()+starts[k+1]);
            r.dedup=secs(t);t=Clock::now();int h=0;U unique=0;
            #pragma omp parallel for reduction(^:h) reduction(+:unique) schedule(dynamic,1)
            for(int k=0;k<K;k++)for(U i=starts[k];i<starts[k+1];i++)if(i==starts[k]||xs[i]!=xs[i-1]){unique++;h^=p3(g,r.c,xs[i]);}
            r.answer=h;r.unique=unique;r.weight=secs(t);
        }else{
            std::sort(xs.begin(),xs.end());xs.erase(std::unique(xs.begin(),xs.end()),xs.end());r.unique=xs.size();r.dedup=secs(t);t=Clock::now();int h=0;
            #pragma omp parallel for reduction(^:h) schedule(dynamic,256)
            for(U i=0;i<xs.size();i++)h^=p3(g,r.c,xs[i]);r.answer=h;r.weight=secs(t);
        }
    }
    if(det && r.visits>r.centers)throw std::runtime_error("deterministic visit budget failed");
    return r;
}
static Result recursive_lv(const Graph&input,U seed){
    Result r;std::mt19937_64 diag(seed^0xB113),rng(seed^0x264264);r.c=diag()&mask(input.n);auto t=Clock::now();
    using Callback=std::function<void(U)>;
    std::function<void(const Graph&,U,const Callback&,bool)>enumerate;
    enumerate=[&](const Graph&g,U c,const Callback&output,bool top){
        int n=g.n;r.nodes++;
        if(n<=8){for(U x=0;x<(U(1)<<n);x++)if(legal(g,c,x)){if(top){r.visits++;r.unique++;}output(x);}return;}
        int left=n/2,right=n-left;
        auto make_aux=[&](int b,int off){Block block{b,off,std::vector<U>(b),{},0,{}};for(int i=0;i<b;i++)for(int j=i+1;j<b;j++)if(rng()&1){block.graph[i]|=U(1)<<j;block.graph[j]|=U(1)<<i;}return block;};
        std::vector<Block>bs{make_aux(left,0),make_aux(right,left)};
        Graph gl{left,bs[0].graph,{}},gr{right,bs[1].graph,{}};columns(gl);columns(gr);
        std::uniform_int_distribution<int> ternary(0,2);
        int shifts[64];for(int i=0;i<n;i++)shifts[i]=ternary(rng);
        enumerate(gl,0,[&](U sl){Center ql=center_for(bs[0],sl);
            enumerate(gr,0,[&](U sr){Center qr=center_for(bs[1],sr),original{ql.z|(qr.z<<left),ql.o|(qr.o<<left)},q;
                for(int i=0;i<n;i++){int old=original.z>>i&1?0:(original.o>>i&1?1:2);int symbol=(old+shifts[i])%3;if(symbol==0)q.z|=U(1)<<i;else if(symbol==1)q.o|=U(1)<<i;}
                if(top)r.centers++;
                Basis b;if(!qsolve(g,q,c,b))return;
                b.stream(n,[&](U z){if(top)r.visits++;U x=(z^q.z)&~q.o,y=z&(q.z|q.o),ox=0,oy=0;
                    for(int i=0;i<n;i++){int s=x>>i&1?1:(y>>i&1?2:0);int old=(s-shifts[i]+3)%3;if(old==1)ox|=U(1)<<i;else if(old==2)oy|=U(1)<<i;}
                    if(owned(bs,original,ox,oy)){if(top)r.unique++;output(x);}
                });
            },false);
        },false);
    };
    enumerate(input,r.c,[&](U x){r.answer^=p3(input,r.c,x);},true);r.enumerate=secs(t);return r;
}
// Disjoint affine branching on a product constraint. No KW cover or ownership.
// This is an experimental decision-tree enumerator, not a new 1.5^n theorem.
static Result branch(const Graph&g,U seed){
    Result r;std::mt19937_64 rng(seed^0xB113);r.c=rng()&mask(g.n);Basis b;int n=g.n;
    auto t=Clock::now();
    std::function<void()>dfs=[&](){r.nodes++;int chosen=-1;U a0=0,b0=0;int pivots[128],cnt=0;bool bad=false;
        for(int i=0;i<n;i++){
            U a=b.reduce(U(1)<<i,0,n),v=b.reduce(g.b[i],int(r.c>>i&1),n);
            if(!a || !v || (a^v)==(U(1)<<n))continue;
            if(!(a&mask(n)) || !(v&mask(n)) || a==v){
                U eq=!(a&mask(n))?v:a;int p=b.add(eq&mask(n),int(eq>>n&1),n);
                if(p==-2){bad=true;break;}if(p>=0){pivots[cnt++]=p;i=-1;chosen=-1;}
            }else if(chosen<0){chosen=i;a0=a;b0=v;}
        }
        if(!bad){
            if(chosen<0)b.stream(n,[&](U x){r.visits++;r.unique++;r.answer^=p3(g,r.c,x);});
            else{
                int p=b.add(a0&mask(n),int(a0>>n&1),n);if(p!=-2)dfs();b.undo(p);
                p=b.add(a0&mask(n),int(a0>>n&1)^1,n);if(p!=-2){int q=b.add(b0&mask(n),int(b0>>n&1),n);if(q!=-2)dfs();b.undo(q);}b.undo(p);
            }
        }while(cnt)b.undo(pivots[--cnt]);
    };dfs();r.enumerate=secs(t);return r;
}
static Result run(const Graph&g,U seed,std::string method,std::string layout,int threads){
    auto t=Clock::now();Result r;
    if(method=="dp-naive")r=dp_naive(g);
    else if(method=="dp-bit")r=dp_bit(g,threads);
    else if(method=="bh-lv"||method=="bh-det")r=bh(g,seed,method=="bh-det");
    else if(method=="bh-det-cache")r=bh_cached(g,seed);
    else if(method=="branch")r=branch(g,seed);
    else if(method=="recursive-lv")r=recursive_lv(g,seed);
    else if(method.rfind("kw-",0)==0||method.rfind("greedy-",0)==0)r=kw(g,seed,method,layout);
    else throw std::runtime_error("unknown method");
    r.total=secs(t);return r;
}
static int brute(const Graph&g){std::vector<int>p(g.n-1);std::iota(p.begin(),p.end(),1);int h=0;do{int prev=0;bool ok=true;for(int x:p){if(!(g.b[prev]>>x&1)){ok=false;break;}prev=x;}if(ok&&(g.b[prev]&1))h^=1;}while(std::next_permutation(p.begin(),p.end()));return h;}
static void selftest(){
    omp_set_num_threads(1);U tests=0;std::vector<std::string>methods={"dp-naive","dp-bit","bh-lv","bh-det","bh-det-cache","branch","recursive-lv","kw-naive","kw-roll","kw-stream","kw-det-owner","kw-det-sort","kw-lv-owner","kw-lv-sort","greedy-det-owner","greedy-det-sort","greedy-lv-sort","greedy-det-lowmem","kw-lv-bucket","kw-det-bucket","kw-sparse-det-bucket"};
    for(int n=2;n<=4;n++)for(U bits=0;bits<(U(1)<<(n*(n-1)));bits++){
        Graph g{n,std::vector<U>(n),{}};int k=0;for(int i=0;i<n;i++)for(int j=0;j<n;j++)if(i!=j){if(bits>>k&1)g.b[i]|=U(1)<<j;k++;}columns(g);int h=brute(g);
        U lastc=0,bhc=0;for(auto method:methods){auto r=run(g,17,method,"geo",1);if(r.answer!=h)throw std::runtime_error("selftest "+method+" n="+std::to_string(n)+" bits="+std::to_string(bits));if(method=="bh-det")bhc=r.c;if(method=="bh-det-cache"&&r.c!=bhc)throw std::runtime_error("BH cached CE");if(method=="kw-naive")lastc=r.c;if((method=="kw-roll"||method=="kw-stream"||method=="kw-det-sort"||method=="kw-sparse-det-bucket")&&r.c!=lastc)throw std::runtime_error("CE mismatch");}tests++;
    }
    for(int n=5;n<=10;n++)for(U seed=1;seed<=12;seed++){
        Graph g=graph(n,seed,seed%3==0?"sparse":"dense");int h=brute(g);
        for(auto method:methods){auto r=run(g,seed,method,"two",1);if(r.answer!=h)throw std::runtime_error("random selftest "+method);if(method.rfind("kw-",0)==0||method=="branch"||method.rfind("bh-",0)==0){U K=0;for(U x=0;x<(U(1)<<n);x++)K+=legal(g,r.c,x);if(K!=r.unique)throw std::runtime_error("P2 count "+method);}}tests++;
    }
    for(int n=11;n<=18;n++)for(U seed=1;seed<=4;seed++){
        auto g=graph(n,seed,"dense");auto h=dp_bit(g,1).answer;auto r=recursive_lv(g,seed);
        if(r.answer!=h)throw std::runtime_error("recursive depth selftest");
        tests++;
    }
    for(int b=1;b<=10;b++){
        auto block=greedy_block(b,0);for(U index=0;index<block.owners.size();index++){
            int owner=block.owners[index];if(owner<0)throw std::runtime_error("greedy uncovered");auto q=block.list[owner];U v=index;
            for(int i=0;i<b;i++){int digit=v%3;v/=3;int prohibited=q.z>>i&1?0:(q.o>>i&1?1:2);if(prohibited==digit)throw std::runtime_error("greedy ownership");}
        }
    }
    omp_set_num_threads(4);
    for(int n=12;n<=20;n+=2)for(U seed=1;seed<=4;seed++){
        auto g=graph(n,seed,"dense");int h=dp_bit(g,1).answer;
        auto bhfull=bh(g,seed,true),bhfast=bh_cached(g,seed);
        if(bhfull.c!=bhfast.c||bhfast.answer!=h)throw std::runtime_error("parallel BH cache selftest");
        U fullc=run(g,seed,"kw-det-sort","two",4).c;
        for(auto method:{"kw-lv-bucket","kw-det-bucket","kw-sparse-det-bucket"}){
            auto r=run(g,seed,method,"two",4);if(r.answer!=h)throw std::runtime_error("parallel bucket selftest");
            if(std::string(method).find("det")!=std::string::npos && r.c!=fullc)throw std::runtime_error("parallel sparse CE");
            U K=0;for(U x=0;x<(U(1)<<n);x++)K+=legal(g,r.c,x);if(K!=r.unique)throw std::runtime_error("parallel bucket count");
        }tests++;
    }
    std::cout<<"{\"selftest\":\"passed\",\"graphs\":"<<tests<<",\"methods\":"<<methods.size()<<"}\n";
}
int main(int argc,char**argv){
    try{
        init_binom();if(argc>1&&std::string(argv[1])=="selftest"){selftest();return 0;}
        if(argc<6){std::cerr<<"bench METHOD N SEED KIND THREADS [two|geo|six]\n";return 2;}
        std::string method=argv[1],kind=argv[4],layout=argc>6?argv[6]:"two";int n=std::stoi(argv[2]),threads=std::stoi(argv[5]);U seed=std::stoull(argv[3]);
        if(n<2||n>40)throw std::runtime_error("n must be 2..40");
        omp_set_num_threads(threads);
        Graph g=graph(n,seed,kind);Result r=run(g,seed,method,layout,threads);rusage usage{};getrusage(RUSAGE_SELF,&usage);
        U hash=1469598103934665603ULL;for(U row:g.b)hash=(hash^row)*1099511628211ULL;
        std::cout<<std::setprecision(9)<<"{\"method\":\""<<method<<"\",\"n\":"<<n<<",\"seed\":"<<seed<<",\"kind\":\""<<kind<<"\",\"threads\":"<<threads<<",\"layout\":\""<<layout<<"\",\"answer\":"<<r.answer<<",\"graph_hash\":"<<hash<<",\"c\":"<<r.c<<",\"centers\":"<<r.centers<<",\"visits\":"<<r.visits<<",\"unique\":"<<r.unique<<",\"nodes\":"<<r.nodes<<",\"table_bytes\":"<<r.tablebytes<<",\"rss_kib\":"<<usage.ru_maxrss<<",\"setup_s\":"<<r.setup<<",\"choose_s\":"<<r.choose<<",\"enumerate_s\":"<<r.enumerate<<",\"dedup_s\":"<<r.dedup<<",\"weight_s\":"<<r.weight<<",\"total_s\":"<<r.total<<"}\n";
    }catch(const std::exception&e){std::cerr<<e.what()<<"\n";return 1;}
}

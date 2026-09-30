// Deterministic, single-thread, polynomial-space implementations.
#define main baseline_main
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wreturn-type"
#include "baseline.cpp"
#pragma GCC diagnostic pop
#undef main
#include "input_families.hpp"

struct PolyResult {
    Result r;
    U local_cache_bytes=0;
    std::string shortcut;
};

template<class F> static bool prefixes_serial(int n,F f){
    for(int k=0;k<=n/2;k++)for(int type=0;type<2;type++){
        int len=n-k,positions=len-type;
        if(positions<k)continue;
        U count=binom[positions][k],s=mask(k);
        for(U j=0;j<count;j++){
            if(!f(Prefix{s|(type?(U(1)<<(len-1)):0),len,k,0}))return false;
            if(k&&j+1<count)s=nextcomb(s);
        }
    }
    return true;
}
static U prefix_count(int n){
    U M=0;for(int k=0;k<=n/2;k++)for(int type=0;type<2;type++)
        if(n-k-type>=k)M+=binom[n-k-type][k];
    return M;
}
// Independent rows give equal conditional contributions. Reconstruct only the
// earlier rows; a dependent current row contributes one signed exact weight.
static U bh_choose_delta(const Graph&g){
    U c=0;
    for(int i=0;i<g.n;i++){
        __int128 delta=0;
        prefixes_serial(g.n,[&](Prefix p){
            if(!(p.x>>i&1))return true;
            Basis b;U active=p.x&mask(i);int fixed=pc(active);
            while(active){int j=low(active);active&=active-1;
                if(b.add(g.b[j]>>p.len,par(g.b[j]&p.x)^int(c>>j&1),p.k)==-2)return true;
            }
            U v=(g.b[i]>>p.len)|(U(par(g.b[i]&p.x))<<p.k);
            while(v&mask(p.k)){int pivot=low(v);if(!b.row[pivot])return true;v^=b.row[pivot];}
            int depth=1+p.k-pc(p.x)+fixed-b.rank;
            __int128 w=__int128(1)<<depth;
            delta+=v?-w:w;
            return true;
        });
        if(delta>0)c|=U(1)<<i;
    }
    return c;
}
static bool bh_enumerate(const Graph&g,U c,U limit,Result&r){
    r.c=c;
    return prefixes_serial(g.n,[&](Prefix p){
        r.centers++;Basis b;U active=p.x;
        while(active){int i=low(active);active&=active-1;
            if(b.add(g.b[i]>>p.len,par(g.b[i]&p.x)^int(c>>i&1),p.k)==-2)return true;
        }
        // An explicit Gray loop permits an immediate abort at the work budget.
        U x=b.particular(p.k),dirs[64];int d=0;
        for(int j=0;j<p.k;j++)if(!b.row[j]){
            U z=U(1)<<j;for(int t=p.k-1;t>=0;t--)if(b.row[t]&&par(b.row[t]&z))z|=U(1)<<t;
            dirs[d++]=z;
        }
        for(U step=0;step<(U(1)<<d);step++){
            if(step)x^=dirs[low(step)];
            if(++r.visits>limit)return false;
            U full=p.x|(x<<p.len);
            if(legal(g,c,full)){r.unique++;r.answer^=p3(g,c,full);}
        }
        return true;
    });
}
static PolyResult bh_poly(const Graph&g){
    PolyResult out;
    auto t=Clock::now();U c=bh_choose_delta(g);out.r.choose=secs(t);t=Clock::now();
    if(!bh_enumerate(g,c,prefix_count(g.n),out.r))throw std::runtime_error("BH CE visit bound");
    out.r.enumerate=secs(t);return out;
}

static U owner_subset(const Block&b,U x,U y){
    U active=x|y;Basis owner;U bits=active;
    while(bits){int i=low(bits);bits&=bits-1;
        U a=(b.graph[i]&active)^(((y>>i)&1)<<i);
        if(owner.add(a,int(y>>i&1),b.b)==-2)throw std::runtime_error("owner inconsistent");
    }
    return owner.particular(b.b);
}
struct PolyCover {
    std::vector<Block> bs;
    U M=1,cache_bytes=0;
    PolyCover(int n,const std::string&layout){
        int remaining=n,off=0;
        while(remaining){
            int b=layout=="wide"?(7*remaining+11)/12:(remaining+1)/2;
            Block block;
            bool reused=false;
            for(const auto&old:bs)if(old.b==b){block=old;block.off=off;reused=true;break;}
            if(!reused){
                block=build_block(b,off,true,0,false);
                if(block.count<=U(n)*n){
                    for(U s=0;s<(U(1)<<b);s++)if(even(block,s))block.list.push_back(center_for(block,s));
                }
                U ternary=1;for(int i=0;i<b&&ternary<=U(n)*n;i++)ternary*=3;
                if(ternary<=U(n)*n){
                    block.owners.resize(ternary);
                    for(U code=0;code<ternary;code++){
                        U v=code,x=0,y=0;
                        for(int i=0;i<b;i++){int s=v%3;v/=3;if(s==1)x|=U(1)<<i;if(s==2)y|=U(1)<<i;}
                        block.owners[code]=int(owner_subset(block,x,y));
                    }
                }
            }
            cache_bytes+=block.list.capacity()*sizeof(Center)+block.owners.capacity()*sizeof(int);
            M*=block.count;bs.push_back(std::move(block));off+=b;remaining-=b;
        }
    }
    template<class F> bool choices(int j,F f)const{
        const auto&b=bs[j];
        if(!b.list.empty()){for(auto q:b.list)if(!f(q))return false;}
        else for(U s=0;s<(U(1)<<b.b);s++)if(even(b,s)&&!f(center_for(b,s)))return false;
        return true;
    }
    bool owns(const Graph&g,U c,Center q,U x)const{
        // Test the cheap, small blocks first. Only compute the needed syndrome
        // coordinates; a rejected candidate never pays for all n coordinates.
        for(auto it=bs.rbegin();it!=bs.rend();++it){
            const auto&b=*it;
            U xx=(x>>b.off)&mask(b.b),yy=0,owner,bits=mask(b.b)^xx;
            while(bits){int i=low(bits);bits&=bits-1;if(par(g.b[b.off+i]&x)^int(c>>(b.off+i)&1))yy|=U(1)<<i;}
            if(!b.owners.empty()){
                U index=0,power=1;
                for(int i=0;i<b.b;i++){index+=power*((xx>>i&1)?1:((yy>>i&1)?2:0));power*=3;}
                owner=U(b.owners[index]);
            }else owner=owner_subset(b,xx,yy);
            if(owner!=((q.z>>b.off)&mask(b.b)))return false;
        }
        return true;
    }
};

static U kw_choose_delta(const Graph&g,const PolyCover&cover){
    U c=0;int n=g.n;
    for(int ell=1;ell<=n;ell++){
        U bas[64]{};int rank=0;__int128 delta=0;U clip=mask(ell),target=U(1)<<(ell-1);
        auto reduce=[&](U v){while(v){int p=low(v);if(!bas[p])break;v^=bas[p];}return v;};
        auto dfs=[&](auto&&self,int j,U rhs)->void{
            // If the decision bit lies in the current column span, every
            // completion contributes equally to both candidates.
            U reduced_target=reduce(target);if(!reduced_target)return;
            if(j==int(cover.bs.size())){
                U residual=reduce(c^rhs);__int128 w=__int128(1)<<(n-rank);
                if(!residual)delta+=w;
                else if(!reduce((c^rhs)^target))delta-=w;
                return;
            }
            const auto&b=cover.bs[j];
            cover.choices(j,[&](Center q){
                int pivots[64],count=0;U next_rhs=rhs;
                for(int k=0;k<b.b;k++){
                    int i=b.off+k;U v;
                    if(q.z>>k&1){v=g.col[i]^(U(1)<<i);next_rhs^=g.col[i]&clip;}
                    else if(q.o>>k&1)v=U(1)<<i;else v=g.col[i];
                    v&=clip;
                    while(v){int p=low(v);if(!bas[p]){bas[p]=v;pivots[count++]=p;rank++;break;}v^=bas[p];}
                }
                self(self,j+1,next_rhs);
                while(count){bas[pivots[--count]]=0;rank--;}
                return true;
            });
        };
        dfs(dfs,0,0);if(delta>0)c|=target;
    }
    return c;
}

// Use x as the parameter. Each block adds local rows to one rollback basis:
// prohibit 0 -> (B_i+e_i)x=1+c_i; prohibit 1 -> x_i=0;
// prohibit 2 -> B_i x=c_i. Inconsistent partial branches are discarded.
static bool kw_enumerate(const Graph&g,const PolyCover&cover,U c,U limit,Result&r){
    Basis basis;bool stopped=false;r.c=c;
    auto dfs=[&](auto&&self,int j,Center global_q)->void{
        if(stopped)return;
        if(j==int(cover.bs.size())){
            r.centers++;
            U x=basis.particular(g.n),dirs[64];int d=0;
            for(int i=0;i<g.n;i++)if(!basis.row[i]){
                U z=U(1)<<i;
                for(int p=g.n-1;p>=0;p--)if(basis.row[p]&&par(basis.row[p]&z))z|=U(1)<<p;
                dirs[d++]=z;
            }
            for(U step=0;step<(U(1)<<d);step++){
                if(step)x^=dirs[low(step)];
                if(++r.visits>limit){stopped=true;return;}
                if(cover.owns(g,c,global_q,x)){r.unique++;r.answer^=p3(g,c,x);}
            }
            return;
        }
        const auto&block=cover.bs[j];
        cover.choices(j,[&](Center q){
            r.nodes++;int pivots[64],count=0;bool consistent=true;
            for(int k=0;k<block.b;k++){
                int i=block.off+k;U a;int rhs;
                if(q.z>>k&1){a=g.b[i]^(U(1)<<i);rhs=1^int(c>>i&1);}
                else if(q.o>>k&1){a=U(1)<<i;rhs=0;}
                else {a=g.b[i];rhs=int(c>>i&1);}
                int p=basis.add(a,rhs,g.n);
                if(p==-2){consistent=false;break;}
                if(p>=0)pivots[count++]=p;
            }
            if(consistent)self(self,j+1,{global_q.z|(q.z<<block.off),global_q.o|(q.o<<block.off)});
            while(count)basis.undo(pivots[--count]);
            return !stopped;
        });
    };
    dfs(dfs,0,{});return !stopped;
}
static PolyResult kw_poly(const Graph&g,const std::string&layout){
    PolyResult out;auto t=Clock::now();PolyCover cover(g.n,layout);
    out.r.setup=secs(t);out.local_cache_bytes=cover.cache_bytes;
    t=Clock::now();U c=kw_choose_delta(g,cover);out.r.choose=secs(t);t=Clock::now();
    if(!kw_enumerate(g,cover,c,cover.M,out.r))throw std::runtime_error("KW CE visit bound");
    out.r.enumerate=secs(t);out.r.centers=cover.M;return out;
}

static PolyResult run_poly(const Graph&g,const std::string&method,const std::string&layout){
    auto t=Clock::now();PolyResult out;auto f=features(g);
    if(!f.strongly_connected)out.shortcut="not_strongly_connected";
    else if(f.bipartite&&f.left.size()!=f.right.size())out.shortcut="unbalanced_bipartite";
    else if(method=="bh-original")out.r=bh(g,0,true);
    else if(method=="kw-original")out.r=kw(g,0,"kw-stream",layout=="half"?"geo":layout);
    else if(method=="bh-ce")out=bh_poly(g);
    else if(method=="kw-ce")out=kw_poly(g,layout);
    else throw std::runtime_error("unknown method");
    out.r.total=secs(t);return out;
}

static void verify_poly(){
    U checked=0;
    auto verify=[&](const Graph&g,int expected){
        for(const auto&layout:{"half","wide"}){
            auto a=bh_poly(g),b=kw_poly(g,layout);
            for(const auto&r:{a,b}){
                if(r.r.answer!=expected)throw std::runtime_error("poly parity");
                U K=0;for(U x=0;x<(U(1)<<g.n);x++)K+=legal(g,r.r.c,x);
                if(K!=r.r.unique)throw std::runtime_error("poly P2");
            }
            if(a.r.c!=bh(g,0,true).c)throw std::runtime_error("BH delta CE");
            if(std::string(layout)=="half"&&b.r.c!=kw(g,0,"kw-stream","geo").c)throw std::runtime_error("KW delta CE");
        }
        checked++;
    };
    for(int n=2;n<=4;n++)for(U bits=0;bits<(U(1)<<(n*(n-1)));bits++){
        Graph g{n,std::vector<U>(n),{}};int k=0;
        for(int i=0;i<n;i++)for(int j=0;j<n;j++)if(i!=j){if(bits>>k&1)g.b[i]|=U(1)<<j;k++;}
        columns(g);verify(g,dp_bit(g,1).answer);
    }
    for(int n=5;n<=12;n++)for(U seed=1;seed<=3;seed++)for(const std::string kind:{"dense","p250","p750","out2","planted050","tournament","bipartite"}){
        auto g=make_input(n,seed,kind);verify(g,dp_bit(g,1).answer);
    }
    for(int n:{8,10,12})for(const std::string kind:{"empty","complete","cycle"}){
        auto g=graph(n,1,kind);verify(g,dp_bit(g,1).answer);
    }
    std::cout<<"{\"status\":\"passed\",\"graphs\":"<<checked<<",\"full_determinization\":true}\n";
}

int main(int argc,char**argv){try{
    init_binom();omp_set_dynamic(0);omp_set_num_threads(1);
    if(argc==2&&std::string(argv[1])=="verify"){verify_poly();return 0;}
    if(argc!=6)throw std::runtime_error("bench_poly METHOD N GRAPH_SEED KIND LAYOUT");
    int n=std::stoi(argv[2]);U seed=std::stoull(argv[3]);std::string method=argv[1],kind=argv[4],layout=argv[5];
    if(n<2||n>52||(layout!="half"&&layout!="wide"))throw std::runtime_error("range/layout");
    Graph g=make_input(n,seed,kind);auto out=run_poly(g,method,layout);const auto&r=out.r;
    U hash=1469598103934665603ULL;for(U row:g.b)hash=(hash^row)*1099511628211ULL;
    rusage usage{};getrusage(RUSAGE_SELF,&usage);
    std::cout<<std::setprecision(10)<<"{\"method\":\""<<method<<"\",\"n\":"<<n<<",\"graph_seed\":"<<seed<<",\"kind\":\""<<kind<<"\",\"layout\":\""<<layout<<"\",\"threads\":1,\"full_determinization\":true,\"graph_hash\":"<<hash<<",\"answer\":"<<r.answer<<",\"c\":"<<r.c<<",\"centers\":"<<r.centers<<",\"visits\":"<<r.visits<<",\"unique\":"<<r.unique<<",\"nodes\":"<<r.nodes<<",\"local_cache_bytes\":"<<out.local_cache_bytes<<",\"shortcut\":\""<<out.shortcut<<"\",\"rss_kib\":"<<usage.ru_maxrss<<",\"setup_s\":"<<r.setup<<",\"choose_s\":"<<r.choose<<",\"enumerate_s\":"<<r.enumerate<<",\"total_s\":"<<r.total<<"}\n";
}catch(const std::exception&e){std::cerr<<e.what()<<"\n";return 1;}}

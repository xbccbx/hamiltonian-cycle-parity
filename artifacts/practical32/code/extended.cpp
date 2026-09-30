// Expanded practical comparison. Baseline is an immutable snapshot of round one.
#define main baseline_main
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wreturn-type"
#include "baseline.cpp"
#pragma GCC diagnostic pop
#undef main

static Graph make_input(int n,U seed,const std::string&kind){
    if(kind=="dense"||kind=="tournament")return graph(n,seed,kind);
    Graph g{n,std::vector<U>(n),{}};std::mt19937_64 rng(seed);
    bool planted=kind.rfind("planted",0)==0;
    if(kind=="out2"){
        for(int i=0;i<n;i++){
            g.b[i]|=U(1)<<((i+1)%n);
            if(n>2){int j;do{j=int(rng()%n);}while(j==i||j==(i+1)%n);g.b[i]|=U(1)<<j;}
        }
    }else{
        int threshold=kind=="bipartite"?500:std::stoi(kind.substr(planted?7:1));
        if(threshold<0||threshold>1000)throw std::runtime_error("density must be 0..1000");
        for(int i=0;i<n;i++)for(int j=0;j<n;j++)if(i!=j){U r=rng();bool ok=kind!="bipartite"||((i<n/2)!=(j<n/2));if(ok&&int(r%1000)<threshold)g.b[i]|=U(1)<<j;}
        if(planted)for(int i=0;i<n;i++)g.b[i]|=U(1)<<((i+1)%n);
        if(kind=="bipartite"&&n%2==0){int h=n/2;for(int i=0;i<h;i++){g.b[i]|=U(1)<<(i+h);g.b[i+h]|=U(1)<<((i+1)%h);}}
    }
    columns(g);return g;
}
struct Features{bool strongly_connected=false,bipartite=true;std::vector<int>left,right;};
static Features features(const Graph&g){
    Features f;int n=g.n;
    auto reach=[&](const std::vector<U>&rows){U visited=1,front=1;while(front){int i=low(front);front&=front-1;U add=rows[i]&~visited;visited|=add;front|=add;}return visited==mask(n);};
    f.strongly_connected=reach(g.b)&&reach(g.col);
    std::vector<int>color(n,-1),queue;
    for(int start=0;start<n;start++)if(color[start]<0){queue.clear();queue.push_back(start);color[start]=0;
        for(size_t j=0;j<queue.size();j++){int i=queue[j];U adj=g.b[i]|g.col[i];while(adj){int v=low(adj);adj&=adj-1;if(color[v]<0){color[v]=color[i]^1;queue.push_back(v);}else if(color[v]==color[i])f.bipartite=false;}}
    }
    for(int i=0;i<n;i++)(color[i]==0?f.left:f.right).push_back(i);
    return f;
}
static U expand_bits(U x,const std::vector<int>&vertices){U result=0;while(x){int i=low(x);x&=x-1;result|=U(1)<<vertices[i];}return result;}
// BH13 Algorithm B: specialize to the detected bipartition, retaining original
// vertex labels so both algorithms use exactly the same randomized diagonal.
static Result bh_bipartite(const Graph&g,U algo_seed,const Features&f){
    Result r;std::mt19937_64 rng(algo_seed^0xB113);r.c=rng()&mask(g.n);
    int l=int(f.left.size()),k=int(f.right.size());std::vector<U>coeff(l);
    for(int i=0;i<l;i++)for(int j=0;j<k;j++)if(g.b[f.left[i]]>>f.right[j]&1)coeff[i]|=U(1)<<j;
    struct alignas(64) Local{U visits=0;int h=0;};std::vector<Local>local(omp_get_max_threads());auto start=Clock::now();
    r.centers=U(1)<<l;
    #pragma omp parallel for schedule(dynamic,128)
    for(U fixed=0;fixed<r.centers;fixed++){
        U x=expand_bits(fixed,f.left);Basis base;bool ok=true;
        for(int j=0;j<k;j++){int i=f.right[j];if(par(g.b[i]&x)^int(r.c>>i&1))base.add(U(1)<<j,0,k);}
        U active=fixed;while(active){int j=low(active);active&=active-1;if(base.add(coeff[j],int(r.c>>f.left[j]&1),k)==-2){ok=false;break;}}
        if(ok){auto&cur=local[omp_get_thread_num()];base.stream(k,[&](U y){cur.visits++;cur.h^=p3(g,r.c,x|expand_bits(y,f.right));});}
    }
    for(auto v:local){r.visits+=v.visits;r.answer^=v.h;}r.unique=r.visits;r.enumerate=secs(start);return r;
}
struct Execution{Result result;Features feature;std::string kernel,shortcut;};
static Execution execute(const Graph&g,U algo_seed,const std::string&method){
    auto start=Clock::now();Execution e;e.feature=features(g);
    if(!e.feature.strongly_connected)e.shortcut="not_strongly_connected";
    else if(e.feature.bipartite&&e.feature.left.size()!=e.feature.right.size())e.shortcut="unbalanced_bipartite";
    else if(method=="bh-practical"){
        e.kernel=e.feature.bipartite?"bh-bipartite":"bh-lv";
        e.result=e.feature.bipartite?bh_bipartite(g,algo_seed,e.feature):bh(g,algo_seed,false);
    }else if(method=="kw-practical"){
        e.kernel="kw-lv-bucket";e.result=kw(g,algo_seed,"kw-lv-bucket","two");
    }else throw std::runtime_error("unknown practical method");
    e.result.total=secs(start);return e;
}
static void verify_extended(){
    init_binom();omp_set_num_threads(4);U checked=0;
    for(int n=4;n<=12;n++)for(U gs=1;gs<=4;gs++)for(const std::string kind:{"p100","p250","dense","p750","p900","planted050","out2","tournament","bipartite"}){
        auto g=make_input(n,gs,kind);int expected=dp_bit(g,1).answer;
        U as=gs*101+17;auto a=execute(g,as,"bh-practical"),b=execute(g,as,"kw-practical");
        if(a.result.answer!=expected||b.result.answer!=expected)throw std::runtime_error("extended parity verification");
        if(a.shortcut.empty()&&(a.result.c!=b.result.c||a.result.unique!=b.result.unique))throw std::runtime_error("extended P2 verification");
        checked++;
    }
    // Exercise labels above bit 40 independently through a unique planted cycle.
    for(int n:{41,44,48,50,52}){
        auto g=graph(n,17,"cycle");auto f=features(g);if(!f.strongly_connected)throw std::runtime_error("wide SCC");
        Basis b;U desired=(U(1)<<(n-1))|7;for(int i=0;i<n;i++)if(b.add(U(1)<<i,int(desired>>i&1),n)==-2)throw std::runtime_error("wide basis");
        if(b.particular(n)!=desired||!legal(g,mask(n),mask(n)))throw std::runtime_error("wide masks");
        checked++;
    }
    std::cout<<"{\"status\":\"passed\",\"cases\":"<<checked<<"}\n";
}
int main(int argc,char**argv){try{
    init_binom();if(argc==2&&std::string(argv[1])=="verify"){verify_extended();return 0;}
    if(argc!=7)throw std::runtime_error("extended METHOD N GRAPH_SEED KIND THREADS ALGORITHM_SEED");
    std::string method=argv[1],kind=argv[4];int n=std::stoi(argv[2]),threads=std::stoi(argv[5]);U gs=std::stoull(argv[3]),as=std::stoull(argv[6]);
    if(n<2||n>52||threads<1||threads>32)throw std::runtime_error("require 2<=n<=52 and 1<=threads<=32");
    omp_set_dynamic(0);omp_set_num_threads(threads);Graph g=make_input(n,gs,kind);auto e=execute(g,as,method);auto&r=e.result;
    U hash=1469598103934665603ULL,arcs=0;for(U row:g.b){hash=(hash^row)*1099511628211ULL;arcs+=pc(row);}rusage usage{};getrusage(RUSAGE_SELF,&usage);
    std::cout<<std::setprecision(10)<<"{\"method\":\""<<method<<"\",\"kernel\":\""<<e.kernel<<"\",\"shortcut\":\""<<e.shortcut<<"\",\"n\":"<<n<<",\"graph_seed\":"<<gs<<",\"algorithm_seed\":"<<as<<",\"kind\":\""<<kind<<"\",\"threads\":"<<threads<<",\"graph_hash\":"<<hash<<",\"arcs\":"<<arcs<<",\"strongly_connected\":"<<(e.feature.strongly_connected?"true":"false")<<",\"bipartite\":"<<(e.feature.bipartite?"true":"false")<<",\"answer\":"<<r.answer<<",\"c\":"<<r.c<<",\"centers\":"<<r.centers<<",\"visits\":"<<r.visits<<",\"unique\":"<<r.unique<<",\"rss_kib\":"<<usage.ru_maxrss<<",\"setup_s\":"<<r.setup<<",\"enumerate_s\":"<<r.enumerate<<",\"dedup_s\":"<<r.dedup<<",\"weight_s\":"<<r.weight<<",\"total_s\":"<<r.total<<"}\n";
}catch(const std::exception&e){std::cerr<<e.what()<<"\n";return 1;}}

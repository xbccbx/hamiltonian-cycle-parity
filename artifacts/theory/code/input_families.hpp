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

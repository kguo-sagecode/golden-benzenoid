// Copyright (c) 2026 Krystal Guo, University of Amsterdam
// SPDX-License-Identifier: MIT

#include <algorithm>
#include <array>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <map>
#include <set>
#include <stdexcept>
#include <string>
#include <tuple>
#include <unordered_map>
#include <utility>
#include <vector>
#include <boost/multiprecision/cpp_int.hpp>

using boost::multiprecision::cpp_int;
using Cell = std::pair<int,int>;
using Vertex = std::pair<int,int>;
using Shape = std::vector<Cell>;

static const std::array<Cell,6> DIRS{{
    {1,0},{-1,0},{0,1},{0,-1},{1,-1},{-1,1}
}};
static const std::array<Vertex,6> VOFF{{
    {1,1},{0,2},{-1,1},{-1,-1},{0,-2},{1,-1}
}};

[[noreturn]] static void fail(const std::string& s) { throw std::runtime_error(s); }

static bool lex_nonnegative(const Cell& p) {
    return p.first > 0 || (p.first == 0 && p.second >= 0);
}
static int hex_distance(const Cell& p) {
    return std::max({std::abs(p.first), std::abs(p.second), std::abs(p.first+p.second)});
}
static Cell addc(const Cell& a, const Cell& b) { return {a.first+b.first,a.second+b.second}; }

static std::array<Vertex,6> cell_vertices(const Cell& c) {
    std::array<Vertex,6> a{};
    int cx=2*c.first+c.second, cy=3*c.second;
    for (int i=0;i<6;i++) a[i]={cx+VOFF[i].first,cy+VOFF[i].second};
    return a;
}
static std::pair<Vertex,Vertex> edge_key(Vertex a, Vertex b) {
    if (b<a) std::swap(a,b);
    return {a,b};
}

struct Mask {
    uint64_t lo=0, hi=0;
    bool test(int i) const { return i<64 ? ((lo>>i)&1ULL) : ((hi>>(i-64))&1ULL); }
    void set(int i) { if(i<64) lo|=1ULL<<i; else hi|=1ULL<<(i-64); }
    bool disjoint(const Mask& o) const { return ((lo&o.lo)==0)&&((hi&o.hi)==0); }
    bool subset_of(const Mask& o) const { return ((lo&~o.lo)==0)&&((hi&~o.hi)==0); }
    bool operator==(const Mask& o) const { return lo==o.lo && hi==o.hi; }
};

static Shape normalized(const Shape& s) {
    int mq=s[0].first, mr=s[0].second;
    for(auto p:s){mq=std::min(mq,p.first);mr=std::min(mr,p.second);}    
    Shape t; t.reserve(s.size());
    for(auto p:s)t.push_back({p.first-mq,p.second-mr});
    std::sort(t.begin(),t.end());
    return t;
}
static Cell rot(Cell p){return {-p.second,p.first+p.second};}
static Cell refl(Cell p){return {p.second,p.first};}
static Shape canonical(const Shape& s){
    Shape cur=s,best; bool first=true;
    for(int k=0;k<6;k++){
        Shape a=normalized(cur), r; r.reserve(cur.size());
        for(auto p:cur)r.push_back(refl(p));
        r=normalized(r);
        if(first||a<best){best=a;first=false;}
        if(r<best)best=r;
        for(auto& p:cur)p=rot(p);
    }
    return best;
}

static std::vector<Shape> rooted_cores_7(){
    std::set<Shape> level; level.insert(Shape{{0,0}});
    for(int n=2;n<=7;n++){
        std::set<Shape> next;
        for(const auto& s:level){
            std::set<Cell> S(s.begin(),s.end()),bd;
            for(auto p:s)for(auto d:DIRS){Cell q=addc(p,d);if(!S.count(q)&&lex_nonnegative(q))bd.insert(q);}
            for(auto q:bd){Shape t=s;t.push_back(q);std::sort(t.begin(),t.end());next.insert(t);}
        }
        level.swap(next);
    }
    return std::vector<Shape>(level.begin(),level.end());
}

struct SupportEntry { Vertex v; cpp_int y; };
struct Pattern { Mask one,zero; uint32_t core=0; std::vector<SupportEntry> support; };
struct LitChild { uint8_t var=0,val=0; uint32_t child=0; };
struct Node { uint32_t pid=0; std::vector<LitChild> a; };

static uint32_t read_u32(std::ifstream& f){
    unsigned char b[4]; f.read(reinterpret_cast<char*>(b),4); if(!f)fail("truncated cover file");
    return uint32_t(b[0])|(uint32_t(b[1])<<8)|(uint32_t(b[2])<<16)|(uint32_t(b[3])<<24);
}
static uint8_t read_u8(std::ifstream& f){char c;f.read(&c,1);if(!f)fail("truncated cover file");return uint8_t((unsigned char)c);}

struct Atlas {
    std::vector<Cell> W7;
    std::map<Cell,int> wi;
    std::vector<Shape> cores;
    std::vector<Pattern> patterns;
};

static Atlas read_atlas(const std::string& path){
    std::ifstream f(path); if(!f)fail("cannot open atlas");
    Atlas A; std::string tok; f>>tok; if(tok!="SEVEN_HEXAGON_OBSTRUCTIONS_V1")fail("bad atlas magic");
    int n;f>>tok>>n;if(tok!="W7"||n!=85)fail("bad W7 header");
    A.W7.resize(n);
    for(int j=0;j<n;j++){int i,q,r;f>>tok>>i>>q>>r;if(tok!="W"||i!=j)fail("bad W line");A.W7[j]={q,r};A.wi[{q,r}]=j;}
    f>>tok>>n;if(tok!="CORES"||n!=3652)fail("bad core header");A.cores.resize(n);
    for(int j=0;j<n;j++){
        int i;f>>tok>>i;if(tok!="C"||i!=j)fail("bad C line");
        Shape s(7);for(auto& p:s)f>>p.first>>p.second;A.cores[j]=s;
    }
    f>>tok>>n;if(tok!="PATTERNS")fail("bad pattern header");A.patterns.resize(n);
    for(int j=0;j<n;j++){
        int pid,di,k;f>>tok>>pid>>di>>k;if(tok!="P"||pid!=j)fail("bad P line");
        Pattern p;p.core=di;
        for(int z=0;z<k;z++){int b;f>>b;if(b<0||b>=85)fail("bad one index");p.one.set(b);}
        f>>k;for(int z=0;z<k;z++){int b;f>>b;if(b<0||b>=85)fail("bad zero index");p.zero.set(b);}
        f>>k;p.support.resize(k);
        for(auto& e:p.support){std::string ys;f>>e.v.first>>e.v.second>>ys;e.y=cpp_int(ys);}
        A.patterns[j]=std::move(p);
    }
    if(!f)fail("atlas parse failure");
    return A;
}

struct CoverData { std::vector<Node> nodes; std::vector<uint32_t> roots; };
static CoverData read_cover(const std::string& path, size_t expected_patterns){
    std::ifstream f(path,std::ios::binary);if(!f)fail("cannot open cover");char magic[8];f.read(magic,8);
    if(std::string(magic,8)!="7HEXCVR1")fail("bad cover magic");
    uint32_t nn=read_u32(f),nr=read_u32(f),np=read_u32(f);
    if(np!=expected_patterns||nr!=3652)fail("cover header count mismatch");
    CoverData C;C.nodes.resize(nn);
    for(uint32_t i=0;i<nn;i++){
        Node x;x.pid=read_u32(f);uint32_t k=read_u32(f);x.a.resize(k);
        for(auto& e:x.a){e.var=read_u8(f);e.val=read_u8(f);e.child=read_u32(f);}
        C.nodes[i]=std::move(x);
    }
    C.roots.resize(nr);for(auto& r:C.roots)r=read_u32(f);
    char extra;if(f.read(&extra,1))fail("extra bytes in cover file");
    return C;
}

static cpp_int sq(const cpp_int& x){return x*x;}

static void verify_pattern(const Atlas& A, const Pattern& p, size_t pid){
    if(p.core>=A.cores.size())fail("pattern core out of range");
    if(!p.one.disjoint(p.zero))fail("pattern has conflicting literals");
    const Shape& C=A.cores[p.core];std::set<Cell> CS(C.begin(),C.end()), collar;
    for(auto c:C)for(auto d:DIRS){Cell z=addc(c,d);if(!CS.count(z)&&lex_nonnegative(z))collar.insert(z);}
    for(auto c:C){auto it=A.wi.find(c);if(it==A.wi.end()||!p.one.test(it->second))fail("pattern omits a core cell");}
    for(int i=0;i<85;i++)if(p.one.test(i)||p.zero.test(i)){
        Cell c=A.W7[i];if(!CS.count(c)&&!collar.count(c))fail("pattern literal outside one-ring");
        if(CS.count(c)&&p.zero.test(i))fail("core cell required absent");
    }
    std::set<Vertex> coreverts;for(auto c:C){auto v=cell_vertices(c);coreverts.insert(v.begin(),v.end());}
    if(p.support.empty())fail("zero witness");
    std::map<Vertex,cpp_int> ymap;int colour=-1;
    for(auto e:p.support){
        if(!coreverts.count(e.v))fail("witness vertex outside core");
        int m=((e.v.second%3)+3)%3;if(m!=1&&m!=2)fail("bad honeycomb colour residue");
        if(colour<0)colour=m;else if(colour!=m)fail("witness uses both colours");
        if(ymap.count(e.v))fail("duplicate witness vertex");ymap[e.v]=e.y;
    }
    cpp_int norm=0;for(auto kv:ymap)norm+=sq(kv.second);if(norm==0)fail("zero witness norm");

    struct EC { bool core=false; std::set<int> covers; };
    std::map<std::pair<Vertex,Vertex>,EC> emap;
    auto addcell=[&](Cell c,bool iscore){
        auto vv=cell_vertices(c);
        int gi=-1;auto it=A.wi.find(c);if(it!=A.wi.end())gi=it->second;
        for(int j=0;j<6;j++){
            auto& e=emap[edge_key(vv[j],vv[(j+1)%6])];
            if(iscore)e.core=true;else {if(gi<0)fail("collar outside W7");e.covers.insert(gi);}        
        }
    };
    for(auto c:C)addcell(c,true);for(auto c:collar)addcell(c,false);

    struct Group { cpp_int fixed=0; std::vector<cpp_int> unknown; };
    std::map<Vertex,Group> groups;
    for(auto kv:ymap){
        Vertex v=kv.first;const cpp_int& coef=kv.second;
        for(auto const& ee:emap){
            Vertex w;bool incident=false;
            if(ee.first.first==v){w=ee.first.second;incident=true;}
            else if(ee.first.second==v){w=ee.first.first;incident=true;}
            if(!incident)continue;
            bool fixed=ee.second.core, unknown=false;
            for(int g:ee.second.covers){if(p.one.test(g))fixed=true;else if(!p.zero.test(g))unknown=true;}
            if(fixed)groups[w].fixed+=coef;else if(unknown)groups[w].unknown.push_back(coef);
        }
    }
    cpp_int total=0;
    for(auto& kv:groups){
        auto& g=kv.second;size_t k=g.unknown.size();if(k>10)fail("impossible unknown degree");
        cpp_int best=0;
        for(size_t m=0;m<(size_t(1)<<k);m++){
            cpp_int s=g.fixed;for(size_t j=0;j<k;j++)if((m>>j)&1)s+=g.unknown[j];
            cpp_int z=s*s;if(z>best)best=z;
        }
        total+=best;
    }
    cpp_int qk=total-2*norm;
    // Exact integer test: with N = y^T y and M = -qhat, the condition
    // phi*N < M is equivalent to 2M-N > 0 and 5N^2 < (2M-N)^2.
    cpp_int N=norm, M=-qk, D=2*M-N;
    if(!(D>0 && 5*N*N<D*D))fail("integer obstruction certificate failed at pattern "+std::to_string(pid));
}

struct StateKey {uint32_t id;Mask o,z;};

static void verify_cover(const Atlas& A, const CoverData& C){
    if(C.nodes.empty())fail("empty cover");
    std::vector<char> seen(C.nodes.size(),0);std::vector<Mask> so(C.nodes.size()),sz(C.nodes.size());
    std::function<void(uint32_t,Mask,Mask)> rec = [&](uint32_t id,Mask one,Mask zero){
        if(id>=C.nodes.size())fail("cover node out of range");const Node& n=C.nodes[id];
        if(n.pid>=A.patterns.size())fail("pattern id out of range");const Pattern& p=A.patterns[n.pid];
        if(!p.one.disjoint(zero)||!p.zero.disjoint(one))fail("node pattern incompatible with cube");
        std::vector<std::pair<uint8_t,uint8_t>> exp;
        for(int b=0;b<85;b++){
            if(p.one.test(b)&&!one.test(b))exp.push_back({uint8_t(b),1});
            else if(p.zero.test(b)&&!zero.test(b))exp.push_back({uint8_t(b),0});
        }
        std::stable_sort(exp.begin(),exp.end(),[](auto a,auto b){
            if(a.second!=b.second)return a.second>b.second;return a.first<b.first;
        });
        if(exp.size()!=n.a.size())fail("node arity does not equal unmatched pattern literals");
        for(size_t j=0;j<exp.size();j++)if(n.a[j].var!=exp[j].first||n.a[j].val!=exp[j].second)fail("node literal order mismatch");
        if(!n.a.empty()){
            if(seen[id]){if(!(so[id]==one&&sz[id]==zero))fail("internal DAG node reached with two cubes");return;}
            seen[id]=1;so[id]=one;sz[id]=zero;
        }
        Mask o=one,z=zero;
        for(auto e:n.a){
            if(e.val){Mask zz=z;zz.set(e.var);rec(e.child,o,zz);o.set(e.var);}else{Mask oo=o;oo.set(e.var);rec(e.child,oo,z);z.set(e.var);}
        }
    };
    for(size_t i=0;i<A.cores.size();i++){
        Mask one,zero;for(auto c:A.cores[i])one.set(A.wi.at(c));rec(C.roots[i],one,zero);
    }
}

// Exact arithmetic in Q(phi), phi^2=phi+1.
static cpp_int iabs(cpp_int x){return x<0?-x:x;}
static cpp_int igcd(cpp_int a,cpp_int b){a=iabs(a);b=iabs(b);while(b!=0){cpp_int r=a%b;a=b;b=r;}return a;}
struct F {
    cpp_int a=0,b=0,d=1; // (a+b phi)/d, d>0
    F()=default;F(long long x):a(x),b(0),d(1){}
    F(cpp_int A,cpp_int B,cpp_int D=1):a(A),b(B),d(D){norm();}
    void norm(){
        if(d==0)fail("zero denominator");if(d<0){a=-a;b=-b;d=-d;}
        if(a==0&&b==0){d=1;return;}cpp_int g=igcd(igcd(a,b),d);if(g!=0&&g!=1){a/=g;b/=g;d/=g;}
    }
};
static F operator+(const F&x,const F&y){return F(x.a*y.d+y.a*x.d,x.b*y.d+y.b*x.d,x.d*y.d);} 
static F operator-(const F&x,const F&y){return F(x.a*y.d-y.a*x.d,x.b*y.d-y.b*x.d,x.d*y.d);} 
static F operator*(const F&x,const F&y){return F(x.a*y.a+x.b*y.b,x.a*y.b+x.b*y.a+x.b*y.b,x.d*y.d);} 
static F inv(const F&x){
    if(x.a==0&&x.b==0)fail("division by zero");cpp_int N=x.a*x.a+x.a*x.b-x.b*x.b;
    return F(x.d*(x.a+x.b),-x.d*x.b,N);
}
static F operator/(const F&x,const F&y){return x*inv(y);} 
static int sign_pair(const cpp_int&a,const cpp_int&b){
    if(a==0&&b==0)return 0;cpp_int u=2*a+b,v=b;
    if(v==0)return u>0?1:-1;if(u==0)return v>0?1:-1;
    if(u>0&&v>0)return 1;if(u<0&&v<0)return -1;
    cpp_int U=u*u,V=5*v*v;
    if(U==V)fail("irrational sign comparison unexpectedly equal");
    if(u>0&&v<0)return U>V?1:-1;
    return V>U?1:-1;
}
static int sgn(const F&x){return sign_pair(x.a,x.b);}

static bool exact_psd(std::vector<std::vector<F>> M){
    while(!M.empty()){
        int n=M.size(),p=-1;
        for(int i=0;i<n;i++){int s=sgn(M[i][i]);if(s<0)return false;if(s>0&&p<0)p=i;}
        if(p<0){for(int i=0;i<n;i++)for(int j=0;j<n;j++)if(sgn(M[i][j])!=0)return false;return true;}
        if(p!=0){std::swap(M[p],M[0]);for(int i=0;i<n;i++)std::swap(M[i][p],M[i][0]);}
        F piv=M[0][0];std::vector<std::vector<F>> N(n-1,std::vector<F>(n-1));
        for(int i=1;i<n;i++)for(int j=1;j<n;j++)N[i-1][j-1]=M[i][j]-M[i][0]*M[0][j]/piv;
        M.swap(N);
    }
    return true;
}

static std::map<Vertex,std::set<Vertex>> carbon_graph(const Shape& s){
    std::set<std::pair<Vertex,Vertex>> E;std::set<Vertex> V;
    for(auto c:s){auto a=cell_vertices(c);V.insert(a.begin(),a.end());for(int i=0;i<6;i++)E.insert(edge_key(a[i],a[(i+1)%6]));}
    std::map<Vertex,std::set<Vertex>> adj;for(auto v:V)adj[v];for(auto e:E){adj[e.first].insert(e.second);adj[e.second].insert(e.first);}return adj;
}
static bool graph_gap_psd(const Shape& s){
    auto adj=carbon_graph(s);std::vector<Vertex> side[2];
    for(auto kv:adj){int m=((kv.first.second%3)+3)%3;if(m==1)side[0].push_back(kv.first);else if(m==2)side[1].push_back(kv.first);else fail("bad carbon coordinate");}
    for(int c=0;c<2;c++){
        int n=side[c].size();std::vector<std::vector<F>> Q(n,std::vector<F>(n));
        for(int i=0;i<n;i++)for(int j=0;j<n;j++){
            int z=0;for(auto w:adj[side[c][i]])if(adj[side[c][j]].count(w))z++;
            Q[i][j]=F(z-(i==j?2:0),i==j?1:0);
        }
        if(!exact_psd(Q))return false;
    }
    return true;
}

static void verify_small_classification(){
    std::set<Shape> level;level.insert(Shape{{0,0}});std::vector<int> expected{0,1,1,3,7,22,82};
    std::vector<std::pair<int,Shape>> survivors;
    for(int h=1;h<=6;h++){
        if((int)level.size()!=expected[h])fail("small polyhex count mismatch at h="+std::to_string(h));
        for(auto const& s:level)if(graph_gap_psd(s))survivors.push_back({h,s});
        if(h==6)break;
        std::set<Shape> next;
        for(auto const& s:level){std::set<Cell>S(s.begin(),s.end()),bd;for(auto p:s)for(auto d:DIRS){Cell q=addc(p,d);if(!S.count(q))bd.insert(q);}for(auto q:bd){Shape t=s;t.push_back(q);next.insert(canonical(t));}}
        level.swap(next);
    }
    if(survivors.size()!=3||survivors[0].first!=1||survivors[1].first!=2||survivors[2].first!=4)fail("small classification did not yield the expected three sizes");
    std::cout<<"SMALL_POLYHEX_COUNTS 1 1 3 7 22 82\n";
    for(auto const& x:survivors){std::cout<<"SURVIVOR h="<<x.first<<" cells=";for(auto p:x.second)std::cout<<"("<<p.first<<","<<p.second<<")";std::cout<<"\n";}
    std::cout<<"SMALL_BENZENOID_CLASSIFICATION_EXACTLY_VERIFIED\n";
}

int main(int argc,char**argv){
    try{
        if(argc!=3){std::cerr<<"usage: "<<argv[0]<<" seven_hexagon_obstructions.txt seven_hexagon_cover.bin\n";return 2;}
        Atlas A=read_atlas(argv[1]);
        std::vector<Cell> W;
        for(int q=-7;q<=7;q++)for(int r=-7;r<=7;r++){Cell p{q,r};if(hex_distance(p)<=7&&lex_nonnegative(p))W.push_back(p);}std::sort(W.begin(),W.end());
        if(W.size()!=85||W!=A.W7)fail("W7 atlas mismatch");
        int w6=0;for(auto p:W)if(hex_distance(p)<=6)w6++;if(w6!=64)fail("W6 count mismatch");
        auto cores=rooted_cores_7();if(cores.size()!=3652||cores!=A.cores)fail("rooted core atlas is not exhaustive");
        for(size_t i=0;i<A.patterns.size();i++)verify_pattern(A,A.patterns[i],i);
        std::cout<<"OBSTRUCTION_PATTERNS_EXACTLY_VERIFIED "<<A.patterns.size()<<"\n";
        CoverData C=read_cover(argv[2],A.patterns.size());verify_cover(A,C);
        std::cout<<"SEVEN_HEXAGON_COVER_EXACTLY_VERIFIED nodes="<<C.nodes.size()<<" roots="<<C.roots.size()<<"\n";
        verify_small_classification();
        std::cout<<"SEVEN_HEXAGON_THEOREM_EXACTLY_VERIFIED\n";
        return 0;
    }catch(const std::exception&e){std::cerr<<"VERIFICATION_FAILED: "<<e.what()<<"\n";return 1;}
}

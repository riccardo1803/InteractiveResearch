#include <bits/stdc++.h>
using namespace std;
static constexpr int MAXN=384;
static constexpr double PI=3.14159265358979323846;
static mt19937_64 rng(2026111601ULL);
static uniform_real_distribution<double> U(0.0,1.0);
struct Pt{double x,y;};
static double torus(Pt a,Pt b){double x=abs(a.x-b.x),y=abs(a.y-b.y);x=min(x,1-x);y=min(y,1-y);return x*x+y*y;}
int main(int argc,char**argv){int M=argc>1?atoi(argv[1]):160;vector<int> ns{50,80,100,150,200,250,300,350}, mus{3,6,9,12,15,18,21,24},ks{2,3,4,5,7,10,14,18,24,30};
cout<<"{\"metadata\":{\"model\":\"Fixed-n 2D hard RGG on unit torus; uniformly planted k-vertex clique; exact vertex-set recovery\",\"seed\":2026111601,\"trialsPerCell\":"<<M<<",\"rFormula\":\"sqrt(mu/(pi*n))\",\"methods\":{\"VD\":\"top k nodes by observed degree, deterministic vertex-ID tie breaking\",\"CN\":\"scan observed edges in deterministic natural-first order; accept first edge with exactly k-2 shared neighbors forming a clique\"}},\"n\":[";
for(int i=0;i<(int)ns.size();i++)cout<<(i?",":"")<<ns[i];cout<<"],\"mu\":[";for(int i=0;i<(int)mus.size();i++)cout<<(i?",":"")<<mus[i];cout<<"],\"k\":[";for(int i=0;i<(int)ks.size();i++)cout<<(i?",":"")<<ks[i];cout<<"],\"data\":[";
for(int ni=0;ni<(int)ns.size();ni++){if(ni)cout<<",";cout<<"[";int n=ns[ni];
 for(int mi=0;mi<(int)mus.size();mi++){if(mi)cout<<",";int mu=mus[mi];vector<int> vd(ks.size(),0),cn(ks.size(),0);
  for(int rep=0;rep<M;rep++){
   vector<Pt> points(n);for(auto& p:points)p={U(rng),U(rng)};
   vector<int> perm(n);iota(perm.begin(),perm.end(),0);shuffle(perm.begin(),perm.end(),rng);
   double r2=(double)mu/(PI*n);int nb=max(1,(int)floor(1.0/sqrt(r2)));
   vector<vector<int>> cells(nb*nb);
   auto cellindex=[&](Pt p){return min(nb-1,int(p.x*nb))*nb+min(nb-1,int(p.y*nb));};
   for(int i=0;i<n;i++)cells[cellindex(points[i])].push_back(i);
   vector<bitset<MAXN>> natural(n);
   vector<int> baseDeg(n);
   vector<pair<int,int>> naturalEdges;
   for(int i=0;i<n;i++){
    int cx=min(nb-1,int(points[i].x*nb)),cy=min(nb-1,int(points[i].y*nb));
    set<int> inspectedCells;
    for(int xx=-1;xx<=1;xx++)for(int yy=-1;yy<=1;yy++){
     int qx=(cx+xx+nb)%nb,qy=(cy+yy+nb)%nb;
     if(!inspectedCells.insert(qx*nb+qy).second)continue;
     for(int j:cells[qx*nb+qy])if(j>i&&torus(points[i],points[j])<r2){natural[i].set(j);natural[j].set(i);baseDeg[i]++;baseDeg[j]++;naturalEdges.emplace_back(i,j);}
    }
   }
   sort(naturalEdges.begin(),naturalEdges.end());
   for(int ki=0;ki<(int)ks.size();ki++){
    int k=ks[ki];vector<bitset<MAXN>> adj=natural;vector<int> deg=baseDeg;
    vector<pair<int,int>> edges=naturalEdges;
    for(int i=0;i<k;i++)for(int j=i+1;j<k;j++){
      int a=perm[i],b=perm[j];if(!adj[a].test(b)){adj[a].set(b);adj[b].set(a);deg[a]++;deg[b]++;edges.emplace_back(a,b);}
    }
    vector<int> order(n);iota(order.begin(),order.end(),0);
    partial_sort(order.begin(),order.begin()+k,order.end(),[&](int a,int b){return deg[a]!=deg[b]?deg[a]>deg[b]:a<b;});
    bitset<MAXN> planted;for(int i=0;i<k;i++)planted.set(perm[i]);
    bool vdok=true;for(int i=0;i<k;i++)if(!planted.test(order[i])){vdok=false;break;}if(vdok)vd[ki]++;
    bool cnok=false;
    for(auto [a,b]:edges){auto shared=(adj[a]&adj[b]);if((int)shared.count()!=k-2)continue;
      bool complete=true;
      for(int x=shared._Find_first();x<n&&complete;x=shared._Find_next(x)){
        auto others=shared;others.reset(x);
        if((others & ~adj[x]).any()){complete=false;break;}
      }
      if(!complete)continue;
      auto candidate=shared;candidate.set(a);candidate.set(b);cnok=(candidate==planted);break;
    }
    if(cnok)cn[ki]++;
   }
  }
  cout<<"[";for(int ki=0;ki<(int)ks.size();ki++){if(ki)cout<<",";cout<<"["<<vd[ki]<<","<<cn[ki]<<"]";}cout<<"]";
 }
 cout<<"]";
}
cout<<"]}";
}

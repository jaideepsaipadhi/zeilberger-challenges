/* Chomp full-box retrograde solver.
   Position = partition (l_1>=...>=l_a), 0<=l_i<=b, must contain (1,1) unless empty.
   Rank: d_k = l_k + (a-k) strictly decreasing; rank = sum C(d_k, a-k+1)  (colex CNS).
   Componentwise-smaller => smaller rank, so a single increasing-rank sweep is a valid DP order.
   Moves from l: pick row i (l_i>=1), threshold t=j-1 (i==1: t>=1; else t>=0), t<=l_i-1:
     mu_k = min(l_k, t) for k>=i.  rank(mu) = P[i-1] + [C(t+m2,t)-C(t+m1-1,t)] + SS[K+1],
     m2=a-i+1, m1=a-K+1, K = max{k: l_k > t}  (hockey-stick collapse of capped rows).
   tbl[rank] = number of winning moves, capped 255; P-position <=> 0. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
typedef unsigned long long u64;
static u64 *Ctab; static int KMAX;
#define C(n,k) Ctab[(u64)(n)*(KMAX+1)+(k)]
static void binit(int n,int kmax){KMAX=kmax;Ctab=calloc((u64)n*(kmax+1),8);
  for(int i=0;i<n;i++){C(i,0)=1;for(int j=1;j<=kmax&&j<=i;j++)C(i,j)=C(i-1,j-1)+((j<=i-1)?C(i-1,j):0);}}
int A,B; u64 NST; uint8_t*tbl;
static inline u64 rank_full(int a2,int b2){ /* bar a2 x b2 inside A x B box */
  u64 r=0; for(int k=1;k<=A;k++){int l=(k<=a2)?b2:0; r+=C(l+A-k,A-k+1);} return r;}
int main(int argc,char**argv){
  A=atoi(argv[1]); B=atoi(argv[2]);
  binit(A+B+3,A+2);
  NST=C(A+B,A);
  tbl=malloc(NST); if(!tbl){fprintf(stderr,"OOM %llu\n",NST);return 1;}
  int l[64]; u64 P[66],SS[66];
  for(int k=0;k<=A+1;k++){l[k]=0;}
  /* iterate ranks 0..NST-1 by colex successor on the combination {d_a<...<d_1} */
  u64 r=0;
  while(1){
    /* prefix ranks P[i]=sum_{k<=i} C(d_k, A-k+1); suffix SS[k]=sum_{j>=k} */
    P[0]=0; for(int k=1;k<=A;k++) P[k]=P[k-1]+C(l[k]+A-k,A-k+1);
    SS[A+1]=0; for(int k=A;k>=1;k--) SS[k]=SS[k+1]+C(l[k]+A-k,A-k+1);
    int cnt=0;
    for(int i=1;i<=A && l[i]>=1;i++){
      int lo=(i==1)?1:0, K=i;
      for(int t=l[i]-1;t>=lo;t--){
        while(K<A && l[K+1]>t) K++;
        int m2=A-i+1, m1=A-K+1;
        u64 rm = P[i-1] + C(t+m2,m2) - C(t+m1-1,m1-1) + SS[K+1];
        if(tbl[rm]==0) cnt++;
      }
    }
    tbl[r]=(cnt>255)?255:cnt;
    /* colex successor: increment smallest d that can move */
    /* d_k = l_k + A - k ; increment l at deepest k where l[k] < l[k-1] (or < B for k=1) */
    int k=A; while(k>=1){ int cap=(k==1)?B:l[k-1]; if(l[k]<cap) break; k--; }
    if(k<1) break;
    l[k]++; for(int j=k+1;j<=A;j++) l[j]=0;
    r++;
  }
  /* report: winning-move counts for all full bars a2 x b2 */
  printf("BOX %d x %d, states %llu\n",A,B,NST);
  for(int a2=1;a2<=A;a2++){for(int b2=1;b2<=B;b2++){
    u64 rr=rank_full(a2,b2); int c=tbl[rr];
    if(c>=2||a2*b2==1) printf("BAR %d x %d : %d winning moves\n",a2,b2,c);
  }}
  /* squares uniqueness */
  for(int s=2;s<=(A<B?A:B);s++){u64 rr=rank_full(s,s); printf("SQ %d: %d\n",s,tbl[rr]);}
  u64 pc=0; for(u64 i=0;i<NST;i++) if(tbl[i]==0) pc++;
  printf("P-positions in box: %llu\n",pc);
  if(argc>3){ /* dump winning moves of bar argv[3] x argv[4] */
    int a2=atoi(argv[3]),b2=atoi(argv[4]);
    for(int k=1;k<=A;k++) l[k]=(k<=a2)?b2:0;
    P[0]=0; for(int k=1;k<=A;k++) P[k]=P[k-1]+C(l[k]+A-k,A-k+1);
    SS[A+1]=0; for(int k=A;k>=1;k--) SS[k]=SS[k+1]+C(l[k]+A-k,A-k+1);
    for(int i=1;i<=A && l[i]>=1;i++){
      int lo=(i==1)?1:0,K=i;
      for(int t=l[i]-1;t>=lo;t--){
        while(K<A && l[K+1]>t) K++;
        int m2=A-i+1,m1=A-K+1;
        u64 rm=P[i-1]+C(t+m2,m2)-C(t+m1-1,m1-1)+SS[K+1];
        if(tbl[rm]==0) printf("  winning bite at row %d col %d\n",i,t+1);
      }
    }
  }
  return 0;
}

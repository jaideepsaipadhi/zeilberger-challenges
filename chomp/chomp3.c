/* v3: bit-packed P/N table (8x less RAM) + early-exit N-detection (scan stops at first
   winning move; full counts recovered for BARS ONLY in a final pass). Atomic OR for
   thread-safe bit writes within a layer; layer barrier orders reads. */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#ifdef _OPENMP
#include <omp.h>
#endif
typedef unsigned long long u64;
static u64 *Ctab; static int KMAX;
#define C(n,k) Ctab[(u64)(n)*(KMAX+1)+(k)]
static void binit(int n,int kmax){KMAX=kmax;Ctab=calloc((u64)n*(kmax+1),8);
  for(int i=0;i<n;i++){C(i,0)=1;for(int j=1;j<=kmax&&j<=i;j++)C(i,j)=C(i-1,j-1)+((j<=i-1)?C(i-1,j):0);}}
int A,B; u64 NST,NW; u64 *bits;
static u64 rk(const int*l){u64 r=0;for(int k=1;k<=A;k++)r+=C(l[k]+A-k,A-k+1);return r;}
static inline int isN(u64 r){return (bits[r>>6]>>(r&63))&1;}
static inline void setN(u64 r){__atomic_or_fetch(&bits[r>>6],1ULL<<(r&63),__ATOMIC_RELAXED);}
/* mode 0: decide P/N with early exit, set bit. mode 1: return full winning-move count. */
static int do_state(const int*l,int mode){
  u64 P[66],SS[66];
  P[0]=0; for(int k=1;k<=A;k++) P[k]=P[k-1]+C(l[k]+A-k,A-k+1);
  SS[A+1]=0; for(int k=A;k>=1;k--) SS[k]=SS[k+1]+C(l[k]+A-k,A-k+1);
  int cnt=0;
  for(int i=1;i<=A && l[i]>=1;i++){
    int lo=(i==1)?1:0,K=i;
    for(int t=l[i]-1;t>=lo;t--){
      while(K<A && l[K+1]>t) K++;
      int m2=A-i+1,m1=A-K+1;
      u64 rm=P[i-1]+C(t+m2,m2)-C(t+m1-1,m1-1)+SS[K+1];
      if(!isN(rm)){ cnt++; if(mode==0){ setN(P[A]); return 1; } }
    }
  }
  return cnt;
}
static void gen(int k,int rem,int cap,int*l){
  if(k>A){ if(rem==0) do_state(l,0); return; }
  if((u64)rem>(u64)cap*(A-k+1)) return;
  int lo=(rem+(A-k))/(A-k+1);
  int hi=(rem<cap)?rem:cap;
  for(int v=hi;v>=lo;v--){ l[k]=v; gen(k+1,rem-v,v,l); }
}
int main(int argc,char**argv){
  A=atoi(argv[1]);B=atoi(argv[2]);
  binit(A+B+3,A+2);
  NST=C(A+B,A); NW=(NST+63)/64;
  bits=calloc(NW,8); if(!bits){fprintf(stderr,"OOM\n");return 1;}
  time_t t0=time(0);
  for(int s=1;s<=A*B;s++){
    int v1lo=(s+A-1)/A; int v1hi=(s<B)?s:B;
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic,1)
#endif
    for(int v1=v1lo;v1<=v1hi;v1++){
      int l[64]; l[1]=v1;
      if(A==1){ if(s==v1) do_state(l,0); continue; }
      gen(2,s-v1,v1,l);
    }
  }
  printf("BOX %d x %d, states %llu, solve %llds\n",A,B,NST,(long long)(time(0)-t0));
  int l[64];
  for(int a2=1;a2<=A;a2++)for(int b2=1;b2<=B;b2++){
    for(int k=1;k<=A;k++)l[k]=(k<=a2)?b2:0;
    int c=do_state(l,1);
    if(c>=2||a2*b2==1)printf("BAR %d x %d : %d winning moves\n",a2,b2,c);
  }
  u64 nn=0; for(u64 w=0;w<NW;w++) nn+=__builtin_popcountll(bits[w]);
  printf("P-positions in box: %llu\n",NST-nn);
  if(argc>3 && argv[3][0]=='T'){ /* TB dump: P/N of all two-block shapes */
    for(int p=1;p<A;p++)for(int q=1;q<=A-p;q++){
      for(int w=1;w<=B;w++){
        printf("TB %d %d %d ",p,q,w);
        for(int t=0;t<=w;t++){
          for(int k=1;k<=A;k++) l[k]=(k<=p)?w:((k<=p+q)?t:0);
          putchar(isN(rk(l))?'N':'P');
        }
        putchar('\n');
      }
    }
    return 0;
  }
  if(argc>4){
    int a2=atoi(argv[3]),b2=atoi(argv[4]);
    for(int k=1;k<=A;k++)l[k]=(k<=a2)?b2:0;
    u64 P[66],SS[66];
    P[0]=0;for(int k=1;k<=A;k++)P[k]=P[k-1]+C(l[k]+A-k,A-k+1);
    SS[A+1]=0;for(int k=A;k>=1;k--)SS[k]=SS[k+1]+C(l[k]+A-k,A-k+1);
    for(int i=1;i<=A&&l[i]>=1;i++){
      int lo=(i==1)?1:0,K=i;
      for(int t=l[i]-1;t>=lo;t--){
        while(K<A&&l[K+1]>t)K++;
        int m2=A-i+1,m1=A-K+1;
        u64 rm=P[i-1]+C(t+m2,m2)-C(t+m1-1,m1-1)+SS[K+1];
        if(!isN(rm))printf("  winning bite at row %d col %d\n",i,t+1);
      }
    }
  }
  return 0;
}

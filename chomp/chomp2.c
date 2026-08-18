/* v2: area-layered, OpenMP-parallel. Same math as chomp.c (validated), new iteration order:
   process states by increasing area; within a layer states are independent (children have
   strictly smaller area), parallelized over the first-row value. */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#ifdef _OPENMP
#include <omp.h>
#endif
typedef unsigned long long u64;
static u64 *Ctab; static int KMAX;
#define C(n,k) Ctab[(u64)(n)*(KMAX+1)+(k)]
static void binit(int n,int kmax){KMAX=kmax;Ctab=calloc((u64)n*(kmax+1),8);
  for(int i=0;i<n;i++){C(i,0)=1;for(int j=1;j<=kmax&&j<=i;j++)C(i,j)=C(i-1,j-1)+((j<=i-1)?C(i-1,j):0);}}
int A,B; u64 NST; uint8_t*tbl;
static inline u64 rk(const int*l){u64 r=0;for(int k=1;k<=A;k++)r+=C(l[k]+A-k,A-k+1);return r;}
static void do_state(const int*l){
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
      if(tbl[rm]==0) cnt++;
    }
  }
  tbl[P[A]]=(cnt>255)?255:cnt;
}
static void gen(int k,int rem,int cap,int*l){
  if(k>A){ if(rem==0) do_state(l); return; }
  if((u64)rem>(u64)cap*(A-k+1)) return;
  int lo=(rem+(A-k))/(A-k+1);
  int hi=(rem<cap)?rem:cap;
  for(int v=hi;v>=lo;v--){ l[k]=v; gen(k+1,rem-v,v,l); }
}
int main(int argc,char**argv){
  A=atoi(argv[1]);B=atoi(argv[2]);
  binit(A+B+3,A+2);
  NST=C(A+B,A);
  tbl=malloc(NST); if(!tbl){fprintf(stderr,"OOM\n");return 1;}
  for(int s=0;s<=A*B;s++){
    int v1lo=(s+A-1)/A; if(v1lo<0)v1lo=0;
    int v1hi=(s<B)?s:B;
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic,1)
#endif
    for(int v1=v1lo;v1<=v1hi;v1++){
      int l[64]; l[1]=v1;
      if(A==1){ if(s==v1){do_state(l);} continue; }
      gen(2,s-v1,v1,l);
    }
    if(s==0){int l[64];for(int k=1;k<=A;k++)l[k]=0;do_state(l);} /* empty state (v1 loop covers s=0,v1=0 already; harmless rewrite) */
  }
  printf("BOX %d x %d, states %llu\n",A,B,NST);
  int l[64];
  for(int a2=1;a2<=A;a2++)for(int b2=1;b2<=B;b2++){
    for(int k=1;k<=A;k++)l[k]=(k<=a2)?b2:0;
    int c=tbl[rk(l)];
    if(c>=2||a2*b2==1)printf("BAR %d x %d : %d winning moves\n",a2,b2,c);
  }
  u64 pc=0;for(u64 i=0;i<NST;i++)if(tbl[i]==0)pc++;
  printf("P-positions in box: %llu\n",pc);
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
        if(tbl[rm]==0)printf("  winning bite at row %d col %d\n",i,t+1);
      }
    }
  }
  return 0;
}

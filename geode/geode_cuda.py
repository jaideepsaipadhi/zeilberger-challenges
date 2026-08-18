#!/usr/bin/env python3
"""CUDA (cupy) backend for the Geode diagonal computation. See CUDA_DESIGN.md.

  python3 geode_cuda.py --selftest              # REQUIRED first: GPU vs validated CPU
  python3 geode_cuda.py --D 4 --M 1000          # acceptance test (compare to Rubine)
  python3 geode_cuda.py --D 5 --M 1000          # the challenge

NOTE: written against cupy; it has NOT been executed on a GPU in the authoring environment.
The selftest is the gate.
"""
import argparse, json, os, math, time
import numpy as np
from sympy import nextprime

KERNEL = r'''
extern "C" __global__
void geode(const unsigned int* __restrict__ primes,
           const unsigned int* __restrict__ pinvs,
           const unsigned int* __restrict__ r2s,      // R^2 mod p
           const unsigned int* __restrict__ Cb,       // [nprimes * (Zdeg+1)]
           const unsigned int* __restrict__ invr,     // [nprimes * (Zdeg+1)]
           const unsigned int* __restrict__ wq,       // [nprimes * n]
           unsigned int* __restrict__ partial,        // [nprimes * n]
           int n, int Zdeg, int M, int d, int sd)
{
    long long gid = blockIdx.x * (long long)blockDim.x + threadIdx.x;
    long long total = (long long)n * gridDim.y;      // nodes x primes handled via gridDim.y
    int pi = blockIdx.y;
    if (gid >= n) return;
    int j = (int)gid;

    unsigned int p = primes[pi], pinv = pinvs[pi], r2 = r2s[pi];

    // ---- Montgomery helpers ----
    #define MM(a,b) ({ unsigned long long t=(unsigned long long)(a)*(b); \
                      unsigned int m=(unsigned int)t*pinv; \
                      unsigned long long u=(t+(unsigned long long)m*p)>>32; \
                      (unsigned int)(u>=p?u-p:u); })
    #define TOM(a) MM((a), r2)

    // ---- node value x_j = j+2 ----
    unsigned int xj = (unsigned int)((j + 2) % p);
    unsigned int one_minus_x = (unsigned int)(((long long)1 - xj) % p + p) % p;

    // ---- g_j(z) = (1-x) * prod_{k=1..d} ((1+z)^k - x), degree sd ----
    unsigned int g[16];
    for (int i=0;i<16;i++) g[i]=0;
    g[0] = TOM(one_minus_x);
    int deg = 0;
    unsigned int mx = TOM(one_minus_x);       // f_0 = 1 - x  (Montgomery)
    for (int k=1;k<=d;k++){
        unsigned int f[8];
        // binomials C(k,i)
        unsigned int c=1;
        for (int i=0;i<=k;i++){
            f[i] = TOM((unsigned int)(c % p));
            c = c*(unsigned int)(k-i)/(unsigned int)(i+1);
        }
        f[0] = mx;                            // replace constant term by (1-x)
        unsigned int ng[16];
        for (int i=0;i<16;i++) ng[i]=0;
        for (int a=0;a<=deg;a++){
            if (!g[a]) continue;
            for (int b=0;b<=k;b++){
                unsigned int t = MM(g[a], f[b]);
                unsigned int s = ng[a+b] + t;
                ng[a+b] = (s>=p)? s-p : s;
            }
        }
        for (int i=0;i<16;i++) g[i]=ng[i];
        deg += k;
    }

    // ---- P_0 = g_0^M ; ig0 = g_0^{-1} (Fermat) ----
    unsigned int base=g[0], P0=TOM(1u), e=(unsigned int)M;
    while(e){ if(e&1) P0=MM(P0,base); base=MM(base,base); e>>=1; }
    unsigned int ig0=TOM(1u); base=g[0]; e=p-2;
    while(e){ if(e&1) ig0=MM(ig0,base); base=MM(base,base); e>>=1; }

    // ---- march the recurrence, contracting on the fly ----
    unsigned int Pb[16];
    for (int i=0;i<16;i++) Pb[i]=0;
    Pb[0]=P0;
    unsigned int acc = MM(P0, Cb[pi*(Zdeg+1)+Zdeg]);
    for (int r=0;r<Zdeg;r++){
        unsigned int s_acc=0;
        int smax = (sd < r+1) ? sd : r+1;
        for (int s=1;s<=smax;s++){
            if (!g[s]) continue;
            long long cm = ((long long)M*s - (r+1) + s) % (long long)p;
            if (cm<0) cm += p;
            if (!cm) continue;
            unsigned int t = MM(g[s], TOM((unsigned int)cm));
            t = MM(t, Pb[s-1]);
            unsigned int q = s_acc + t; s_acc = (q>=p)? q-p : q;
        }
        unsigned int Pn = MM(MM(s_acc, ig0), invr[pi*(Zdeg+1)+r]);
        for (int i=sd;i>0;i--) Pb[i]=Pb[i-1];
        Pb[0]=Pn;
        unsigned int t = MM(Pn, Cb[pi*(Zdeg+1)+Zdeg-(r+1)]);
        unsigned int q = acc + t; acc = (q>=p)? q-p : q;
    }
    partial[pi*(long long)n + j] = MM(acc, wq[pi*(long long)n + j]);
}
'''

def host_tables(M,D,p):
    """invr, Cb, wq for one prime (numpy). Mirrors the validated CPU code."""
    d=D-1; sd=d*(d+1)//2; n=D*M+1; Zdeg=D*M
    E=(sum(range(2,D+2))-sd)*M+2
    inv=lambda a: pow(int(a)%p,p-2,p)
    invr=np.array([inv(r+1) for r in range(Zdeg+1)],dtype=np.int64)
    Cb=np.zeros(Zdeg+1,dtype=np.int64); Cb[0]=1
    for s in range(1,Zdeg+1): Cb[s]=Cb[s-1]*((E-s+1)%p)%p*int(invr[s-1])%p
    xs=(np.arange(n,dtype=np.int64)+2)%p
    Pi=np.zeros(n+1,dtype=np.int64); Pi[0]=1
    for deg,xl in enumerate(xs):
        old=Pi.copy(); Pi[:]=0
        Pi[1:deg+2]=old[0:deg+1]
        Pi[0:deg+1]=(Pi[0:deg+1]-old[0:deg+1]*int(xl))%p
    invi=np.array([inv(i+1) for i in range(n)],dtype=np.int64)
    rem=np.full(n,int(Pi[n])%p,dtype=np.int64); integ=(rem*invi[n-1])%p
    for i in range(n-2,-1,-1):
        rem=(int(Pi[i+1])+rem*xs)%p
        integ=(integ+rem*invi[i])%p
    dv=np.ones(n,dtype=np.int64)
    for l in range(n):
        diff=(xs-int(xs[l]))%p; diff[l]=1
        dv=(dv*diff)%p
    wq=(integ*np.array([inv(v) for v in dv],dtype=np.int64))%p
    return invr,Cb,wq,n,Zdeg,sd,d

def gpu_run(M,D,primes):
    import cupy as cp
    d=D-1; sd=d*(d+1)//2; n=D*M+1; Zdeg=D*M
    mod=cp.RawModule(code=KERNEL,options=('-std=c++11',))
    ker=mod.get_function('geode')
    P=len(primes)
    invr=np.zeros((P,Zdeg+1),dtype=np.uint32); Cb=np.zeros((P,Zdeg+1),dtype=np.uint32)
    wq=np.zeros((P,n),dtype=np.uint32); pinvs=np.zeros(P,dtype=np.uint32); r2s=np.zeros(P,dtype=np.uint32)
    for k,p in enumerate(primes):
        a,b,c,_,_,_,_=host_tables(M,D,p)
        # Montgomery domain conversion happens on device; tables uploaded in normal domain
        invr[k]=a%p; Cb[k]=b%p; wq[k]=c%p
        R=1<<32
        pinvs[k]=(-pow(p,-1,R))%R
        r2s[k]=(R*R)%p
    # NOTE: kernel converts tables with TOM() only where needed; Cb/invr/wq are used inside MM(),
    # so upload them already in Montgomery form:
    for k,p in enumerate(primes):
        R=1<<32
        invr[k]=[(int(v)*R)%p for v in invr[k]]
        Cb[k]  =[(int(v)*R)%p for v in Cb[k]]
        wq[k]  =[(int(v)*R)%p for v in wq[k]]
    dprimes=cp.asarray(np.array(primes,dtype=np.uint32)); dpinv=cp.asarray(pinvs); dr2=cp.asarray(r2s)
    dCb=cp.asarray(Cb.ravel()); dinvr=cp.asarray(invr.ravel()); dwq=cp.asarray(wq.ravel())
    part=cp.zeros(P*n,dtype=cp.uint32)
    threads=256; blocks=((n+threads-1)//threads, P)
    ker((blocks[0],blocks[1]),(threads,),(dprimes,dpinv,dr2,dCb,dinvr,dwq,part,
        np.int32(n),np.int32(Zdeg),np.int32(M),np.int32(d),np.int32(sd)))
    cp.cuda.Stream.null.synchronize()
    res=[]
    host=cp.asnumpy(part).reshape(P,n).astype(object)
    for k,p in enumerate(primes):
        R=1<<32; Rinv=pow(R,-1,p)
        s=int(sum(int(v) for v in host[k])%p)
        s=s*Rinv%p                       # leave Montgomery domain
        fd=1
        for i in range(1,D*M+1): fd=fd*i%p
        fM=1
        for i in range(1,M+1): fM=fM*i%p
        res.append(s*fd%p*pow(pow(fM,D,p),p-2,p)%p)
    return res

if __name__=='__main__':
    ap=argparse.ArgumentParser()
    ap.add_argument('--selftest',action='store_true'); ap.add_argument('--D',type=int,default=5)
    ap.add_argument('--M',type=int,default=1000); ap.add_argument('--bits',type=int,default=30)
    a=ap.parse_args()
    if a.selftest:
        from geode_fast import G_modp_vec
        ps=[int(nextprime(2**30)),int(nextprime(2**30+1000))]
        ok=True
        for D in [4,5]:
            for M in [2,3,5,8]:
                cpu=[G_modp_vec(M,D,p) for p in ps]
                gpu=gpu_run(M,D,ps)
                m=all(int(x)==int(y) for x,y in zip(cpu,gpu)); ok=ok and m
                print("  D=%d M=%d: CPU %s  GPU %s  MATCH=%s"%(D,M,cpu,gpu,m))
        print("SELFTEST", "PASS" if ok else "FAIL")
    else:
        lg=lambda k: math.lgamma(k+1)/math.log(10)
        dg=lg(a.M*sum(range(2,a.D+2))+2)-lg(a.M*sum(range(2,a.D+2))+2-a.D*a.M)-lg(a.M+1)-(a.D-1)*lg(a.M)
        need=int(dg*math.log2(10)/a.bits)+6
        print("~%.0f digits -> %d primes"%(dg,need))
        primes=[]; p=int(nextprime(2**a.bits))
        while len(primes)<need: primes.append(p); p=int(nextprime(p))
        t0=time.time(); res=gpu_run(a.M,a.D,primes); print("kernel+host %.1fs"%(time.time()-t0))
        Mtot=1
        for p in primes: Mtot*=p
        acc=0
        for p,r in zip(primes,res):
            Mi=Mtot//p; acc=(acc+r*Mi*pow(Mi,-1,p))%Mtot
        if acc>Mtot//2: acc-=Mtot
        open('G_D%d_M%d_gpu.txt'%(a.D,a.M),'w').write(str(acc)+"\n")
        print("digits:",len(str(abs(acc)))); print("first 60:",str(acc)[:60])

#!/usr/bin/env python3
"""Production driver: exact G(M,...,M) via multi-prime CRT. Resumable.

  python3 run_geode.py --D 4 --M 1000            # ACCEPTANCE TEST (compare to Rubine)
  python3 run_geode.py --D 5 --M 1000            # the $200 challenge
  python3 run_geode.py --D 5 --M 1000 --gpu      # cupy backend if available

State in state_D{D}_M{M}.json (prime -> residue). Safe to kill/restart.
"""
import argparse, json, os, time, math
from sympy import nextprime
from geode_fast import G_modp_vec

def digits_estimate(M,D):
    lg=lambda n: math.lgamma(n+1)/math.log(10)
    d=D-1; base_k=M*sum(range(2,D+2))+2; N=D*M+1
    return lg(base_k)-lg(base_k-N+1)-lg(M+1)-d*lg(M)

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--D',type=int,default=5); ap.add_argument('--M',type=int,default=1000)
    ap.add_argument('--bits',type=int,default=30); ap.add_argument('--gpu',action='store_true')
    ap.add_argument('--extra',type=int,default=5,help='safety primes beyond the bound')
    a=ap.parse_args()
    dg=digits_estimate(a.M,a.D)
    need=int(dg*math.log2(10)/a.bits)+1+a.extra
    print("D=%d M=%d : answer ~%.0f digits -> %d primes of %d bits"%(a.D,a.M,dg,need,a.bits),flush=True)
    sf='state_D%d_M%d.json'%(a.D,a.M)
    st=json.load(open(sf)) if os.path.exists(sf) else {}
    p=int(nextprime(2**a.bits))
    t0=time.time(); done=len(st)
    while len(st)<need:
        if str(p) not in st:
            r=G_modp_vec(a.M,a.D,p)
            st[str(p)]=int(r)
            json.dump(st,open(sf,'w'))
            n=len(st)
            if n%10==0 or n==done+1:
                el=time.time()-t0
                print("  %d/%d primes  %.1fs elapsed  ETA %.1f min"%(n,need,el,(need-n)*el/max(n-done,1)/60),flush=True)
        p=int(nextprime(p))
    # CRT
    print("reconstructing...",flush=True)
    Mtot=1; res=0
    items=[(int(k),v) for k,v in st.items()][:need]
    for pk,_ in items: Mtot*=pk
    for pk,rk in items:
        Mi=Mtot//pk
        res=(res+rk*Mi*pow(Mi,-1,pk))%Mtot
    if res>Mtot//2: res-=Mtot
    out='G_D%d_M%d.txt'%(a.D,a.M)
    open(out,'w').write(str(res)+"\n")
    print("G(%s) has %d digits; written to %s"%(",".join([str(a.M)]*a.D),len(str(abs(res))),out))
    print("first 60 digits:",str(res)[:60])
    print("last 30 digits:",str(res)[-30:])

if __name__=='__main__': main()

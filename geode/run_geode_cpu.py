#!/usr/bin/env python3
"""Production CPU driver: batched primes + multiprocessing + resumable state.

  python3 run_geode_cpu.py --D 4 --M 1000 --workers 8    # acceptance test vs Rubine
  python3 run_geode_cpu.py --D 5 --M 1000 --workers 8    # the $200 challenge

Each worker process handles batches of --batch primes independently (no shared state);
state_D{D}_M{M}.json accumulates prime->residue and is safe to kill/restart."""
import argparse, json, os, time, math, sys
sys.set_int_max_str_digits(0)   # Python 3.11+: lift the 4300-digit int->str cap
# our workload is pure elementwise integer numpy: BLAS threads only cause harm at high worker counts
for v in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS'):
    os.environ.setdefault(v,'1')
from multiprocessing import Pool
from sympy import nextprime

def batch_job(args):
    M,D,ps=args
    from geode_batch import G_modp_batch
    return ps, G_modp_batch(M,D,ps)

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--D',type=int,default=5); ap.add_argument('--M',type=int,default=1000)
    ap.add_argument('--workers',type=int,default=os.cpu_count())
    ap.add_argument('--batch',type=int,default=6); ap.add_argument('--extra',type=int,default=15)
    a=ap.parse_args()
    # rigorous digit bound: |G| <= (a_0/(M!)^d) sum_i i! C(dM,i)/(M+1+i)!  [verified vs exact]
    d=a.D-1; lg=lambda n: math.lgamma(n+1)
    bk=a.M*sum(range(2,a.D+2))+2; bk1=a.M*sum(range(1,a.D+1))+1
    terms=[lg(i)+lg(d*a.M)-lg(i)-lg(d*a.M-i)-lg(a.M+1+i) for i in range(d*a.M+1)]
    mx=max(terms)
    dg=(lg(bk)-lg(1+bk1)+mx+math.log(sum(math.exp(t-mx) for t in terms))-d*lg(a.M))/math.log(10)
    need=int(dg*math.log2(10)/30)+1+a.extra
    print("D=%d M=%d: <=~%.0f digits -> %d primes; workers=%d batch=%d"%(a.D,a.M,dg,need,a.workers,a.batch),flush=True)
    sf='state_D%d_M%d.json'%(a.D,a.M)
    st=json.load(open(sf)) if os.path.exists(sf) else {}
    primes=[]; p=2**30
    while len(primes)<need:
        p=int(nextprime(p))
        if str(p) not in st: primes.append(p)
        # count already-done primes toward 'need'
    primes=primes[:max(0,need-len(st))]
    batches=[(a.M,a.D,primes[i:i+a.batch]) for i in range(0,len(primes),a.batch)]
    t0=time.time(); done=0
    if batches:
        with Pool(a.workers) as pool:
            for ps,rs in pool.imap_unordered(batch_job,batches):
                for q,r in zip(ps,rs): st[str(q)]=int(r)
                json.dump(st,open(sf,'w'))
                done+=len(ps)
                el=time.time()-t0
                print("  %d/%d new primes  %.0fs  ETA %.1f min"%(done,len(primes),el,(len(primes)-done)*el/done/60),flush=True)
    print("reconstructing from %d residues..."%len(st),flush=True)
    items=sorted(((int(k),v) for k,v in st.items()))[:need]
    Mtot=1
    for q,_ in items: Mtot*=q
    res=0
    for q,r in items:
        Mi=Mtot//q
        res=(res+r*Mi*pow(Mi%q,-1,q))%Mtot
    if res>Mtot//2: res-=Mtot
    out='G_D%d_M%d.txt'%(a.D,a.M)
    open(out,'w').write(str(res)+"\n")
    print("G(%s): %d digits -> %s"%(",".join([str(a.M)]*a.D),len(str(abs(res))),out))
    print("first 60:",str(res)[:60]); print("last 30:",str(res)[-30:])

if __name__=='__main__': main()

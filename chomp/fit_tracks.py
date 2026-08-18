"""v2: fit rotation laws; tail-only error windows; predictivity gate; ground-truth backtest;
collision scan uses ACTUAL data where available, predictions only beyond it."""
import sys, glob, math
import numpy as np
files=sys.argv[1:] or sorted(glob.glob("tb_*.out"))
tb={}
for f in files:
    for line in open(f):
        if line.startswith('TB'):
            _,p,q,w,s=line.split(); key=(int(p),int(q),int(w))
            if key not in tb or len(s)>len(tb[key]): tb[key]=s
maxw={}
for (p,q,w) in tb: maxw[(p,q)]=max(maxw.get((p,q),0),w)
tracks={}; fits={}
print("family | cells | slope | window(tail) | spacing | predictive?")
for (p,q),B in sorted(maxw.items()):
    tr=[(w,t) for w in range(3,B+1) for t in range(0,w) if tb.get((p,q,w),'')[t:t+1]=='P']
    tracks[(p,q)]=(tr,B)
    if len(tr)<12: print("(%d,%d) | %4d | too few cells"%(p,q,len(tr))); continue
    ws=np.array([w for w,_ in tr],dtype=float); n=np.arange(1,len(ws)+1,dtype=float)
    h=len(ws)//3
    A=np.polyfit(n[h:],ws[h:],1)
    errt=ws[h:]-np.polyval(A,n[h:])
    lo,hi=errt.min(),errt.max(); width=hi-lo
    pred = width < 0.75*A[0]
    fits[(p,q)]=(A[0],A[1],lo,hi,pred)
    print("(%d,%d) | %4d | %.5f | [%+.2f,%+.2f] w=%.2f | %.2f | %s"%(p,q,len(tr),A[0],lo,hi,width,A[0],"YES" if pred else "no"))
def actual_hit(pq,b):
    p,q=pq
    s=tb.get((p,q,b))
    return None if s is None else ('P' in s[:b])
def pred_hit(pq,b,slack=0.4):
    a0,a1,lo,hi,ok=fits[pq]
    if not ok: return False
    ns=(b-a1)/a0
    for n in (math.floor(ns),math.ceil(ns)):
        if n>=1 and lo-slack<=b-(a0*n+a1)<=hi+slack: return True
    return False
print("\n== BACKTEST (predictive families, covered range: predictions vs truth) ==")
for pq,(a0,a1,lo,hi,ok) in sorted(fits.items()):
    if not ok: continue
    tr,B=tracks[pq]
    tp=fp=fn=0
    for b in range(12,B+1):
        a=actual_hit(pq,b); pr=pred_hit(pq,b)
        if a and pr: tp+=1
        elif pr and not a: fp+=1
        elif a and not pr: fn+=1
    prec=tp/(tp+fp) if tp+fp else 0; rec=tp/(tp+fn) if tp+fn else 0
    print("  %s: precision %.2f recall %.2f (tp %d fp %d fn %d)"%(pq,prec,rec,tp,fp,fn))
print("\n== COLLISIONS ==")
for a in range(6,12):
    fams=[(p,a-p) for p in range(1,a) if (p,a-p) in tracks]
    if len(fams)<3: continue
    Bmin=min(tracks[f][1] for f in fams)
    data_c=[]
    for b in range(6,Bmin+1):
        h=[f for f in fams if actual_hit(f,b)]
        if len(h)>=2: data_c.append((b,len(h),h))
    print("a=%d (data to w=%d): actual multi-hits: %s"%(a,Bmin,
          [(b,c) for b,c,_ in data_c] or "none"))
    trip=[(b,c,h) for b,c,h in data_c if c>=3]
    if trip: print("   *** ACTUAL TRIPLE(S): %s ***"%trip)
    predf=[f for f in fams if f in fits and fits[f][4]]
    if len(predf)>=3:
        cands=[b for b in range(Bmin+1,Bmin+300) if sum(pred_hit(f,b) for f in predf)>=3]
        print("   beyond-data 3-way candidates (predictive fams %s): %s"%(predf,cands[:12] or "none"))
    else:
        print("   only %d predictive families beyond data -- no honest extrapolation yet"%len(predf))

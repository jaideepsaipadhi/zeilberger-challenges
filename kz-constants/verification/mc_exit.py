"""Monte Carlo certification of the B/D targets for the ACTUAL tilted run-MAP:
   (1) P(tau > n) ~ n^{-3/2}  (p=3 for the pi/3 cone)
   (2) conditioned endpoint angle ~ sin(3*phi) after isotropization by (5M)^{-1/2}."""
import numpy as np
rng=np.random.default_rng(7)
v=np.array([[1,0],[-1,1],[0,-1]],dtype=np.int64)
NPATH=2_000_000
# start slightly inside the quadrant to avoid instant death at the apex
start=np.array([2,1])
pos=np.tile(start,(NPATH,1)).astype(np.int64)
state=rng.integers(0,3,NPATH)
alive=np.ones(NPATH,bool)
checkpoints=[8,16,32,64,128,256,512,1024]
surv={}
KMAX=1100
endpts=None
for k in range(1,KMAX+1):
    idx=np.where(alive)[0]
    if len(idx)==0: break
    st=state[idx]
    # next direction: uniform over the other two
    r=rng.integers(0,2,len(idx))
    nxt=(st+1+r)%3
    # run length: R-1 ~ Geometric(1/2) on {1,2,...}
    R=(rng.geometric(0.5,len(idx))+1).astype(np.int64)
    pos[idx]+=R[:,None]*v[nxt]
    state[idx]=nxt
    # kill on leaving quadrant  (tandem cone: x>=0, y>=0)
    dead=(pos[idx,0]<0)|(pos[idx,1]<0)
    alive[idx[dead]]=False
    if k in checkpoints:
        surv[k]=alive.sum()/NPATH
        if k==512: endpts=pos[alive].copy()
ks=np.array(checkpoints[:len(surv)]); ps=np.array([surv[k] for k in checkpoints[:len(surv)]])
# fit slope on the tail
sl=np.polyfit(np.log(ks[3:]),np.log(ps[3:]),1)[0]
print("survival:", {int(k):float(f"{p:.3e}") for k,p in surv.items()})
print("tail slope of log P(tau>k) vs log k: %.4f   (prediction: -1.5)"%sl)
# angular profile at k=512 survivors, isotropized by (5M)^{-1/2}
M=np.array([[2/3,-1/3],[-1/3,2/3]])
S=5*M
w,U=np.linalg.eigh(S)
T=U@np.diag(w**-0.5)@U.T          # T S T = I
z=endpts@T.T
phi=np.arctan2(z[:,1],z[:,0])
# the transformed cone is a wedge; find its actual angular extent from data, then compare shape
lo,hi=np.quantile(phi,[0.001,0.999])
H,edges=np.histogram(phi,bins=24,range=(lo,hi),density=True)
mid=(edges[1:]+edges[:-1])/2
# normalize sin(3*(phi-phi0)) over the wedge [phi0, phi0+pi/3]
phi0=lo
pred=np.sin(np.clip(3*(mid-phi0),0,np.pi))
pred=pred/np.trapezoid(pred,mid)
err=np.max(np.abs(H-pred))/np.max(pred)
print("wedge measured: [%.4f, %.4f], width %.4f  (prediction width pi/3 = %.4f)"%(lo,hi,hi-lo,np.pi/3))
print("angular density vs sin(3(phi-phi0)): max relative deviation %.3f"%err)
for a,b,c in zip(mid[::4],H[::4],pred[::4]):
    print("   phi=%.3f  data=%.3f  sin3=%.3f"%(a,b,c))

import sys
def load(f):
    tb={}
    for line in open(f):
        if line.startswith('TB'):
            _,p,q,w,s=line.split(); tb[(int(p),int(q),int(w))]=s
    return tb
def track(tb,p,q,B):
    return [(w,t) for w in range(1,B+1) for t in range(0,w) if tb.get((p,q,w),'')[t:t+1]=='P']
def balance_test(word,maxlen=None):
    n=len(word); maxlen=maxlen or n//2
    worst=0
    for L in range(1,maxlen+1):
        cs=[sum(word[i:i+L]) for i in range(n-L+1)]
        spread=max(cs)-min(cs)
        worst=max(worst,spread)
        if spread>1: return False,L,spread
    return True,None,worst
def analyze(tb,p,q,B,name):
    tr=track(tb,p,q,B)
    tr=[c for c in tr if c[0]>2]         # drop degenerate small cells
    if len(tr)<8: print("%s: only %d cells, skip"%(name,len(tr))); return
    ws=[w for w,_ in tr]; ts=[t for _,t in tr]
    # single-valued check
    sv = len(set(ws))==len(ws)
    wg=[ws[i+1]-ws[i] for i in range(len(ws)-1)]
    tg=[ts[i+1]-ts[i] for i in range(len(ts)-1)]
    # gap alphabet
    wa=sorted(set(wg)); ta=sorted(set(tg))
    res=""
    if len(wa)<=2:
        word=[0 if g==wa[0] else 1 for g in wg]
        ok,L,sp=balance_test(word)
        res="w-gaps alphabet %s BALANCED=%s"%(wa,ok if ok else "NO(len %d spread %d)"%(L,sp))
    else:
        res="w-gaps alphabet %s (not 2-letter -> not Sturmian as-is)"%wa
    slope=(ws[-1]-ws[0])/(ts[-1]-ts[0]) if ts[-1]!=ts[0] else float('inf')
    print("%s: %d cells, single-valued=%s, slope w/t=%.5f, %s"%(name,len(tr),sv,slope,res))
    print("   w-gaps:", "".join(str(g) for g in wg[:70]))
tb3=load('tb3.out'); tb4=load('tb4.out')
analyze(tb3,2,1,300,"(2,1) [3-row, w<=300]")
analyze(tb3,1,2,300,"(1,2) [3-row, w<=300]")
analyze(tb4,3,1,150,"(3,1) [4-row, w<=150]")
analyze(tb4,2,2,150,"(2,2) [4-row, w<=150]")
analyze(tb4,1,3,150,"(1,3) [4-row, w<=150]")

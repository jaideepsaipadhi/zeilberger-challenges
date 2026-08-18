tb={}
for line in open('tb7.out'):
    if not line.startswith('TB'): continue
    _,p,q,w,s=line.split()
    tb[(int(p),int(q),int(w))]=s
A=7; B=max(w for (_,_,w) in tb)
ok=True
for w in range(1,B+1):
    s=tb[(1,1,w)]
    for t in range(0,w+1):
        expect='P' if (w==t+1 or (w==1 and t==0)) else 'N'
        if s[t]!=expect: ok=False; print("2-row law FAIL",w,t,s[t])
print("2-row classical law:","PASS" if ok else "FAIL")
print("\n-- diagonals t=w-d: eventual period (tail from index 8) --")
for p in range(1,A):
    for q in range(1,A-p+1):
        pats=[]
        for d in range(1,12):
            seq=''.join(tb[(p,q,w)][w-d] for w in range(max(d,1),B+1))
            tail=seq[8:]
            per=0
            for T in range(1,13):
                if len(tail)>=2*T and all(tail[i]==tail[i+T] for i in range(len(tail)-T)): per=T;break
            pats.append("d%d:%s"%(d,per if per else "?"))
        print("(p=%d,q=%d) %s"%(p,q," ".join(pats)))
print("\n-- P-cells (w,t), t<w, per (p,q) with p+q<=7, p,q>=1 --")
for p in range(1,A):
    for q in range(1,A-p+1):
        cells=[(w,t) for w in range(1,B+1) for t in range(0,w) if tb[(p,q,w)][t]=='P']
        if cells: print("(p=%d,q=%d):"%(p,q),cells if len(cells)<28 else str(cells[:24])+" ...(%d total)"%len(cells))

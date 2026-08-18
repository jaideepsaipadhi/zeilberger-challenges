import subprocess,sys
def cvals(A,B):
    out=subprocess.run(["./chomp",str(A),str(B)],capture_output=True,text=True).stdout
    bars={};pc=None
    for line in out.splitlines():
        if line.startswith("BAR"):
            p=line.split(); bars[(int(p[1]),int(p[3]))]=int(p[5])
        if line.startswith("P-positions"): pc=int(line.split()[-1])
    return bars,pc,out
def bvals(A,B):
    out=subprocess.run(["python3","brute.py",str(A),str(B)],capture_output=True,text=True).stdout
    bars={};pc=None
    for line in out.splitlines():
        if line.startswith("BRUTE") and "P-pos" not in line:
            p=line.split(); bars[(int(p[1]),int(p[3]))]=int(p[5])
        if "P-positions" in line: pc=int(line.split()[-1])
    return bars,pc
allok=True
for (A,B) in [(3,3),(4,4),(5,5),(4,6),(2,7),(6,6)]:
    cb,cp,cout=cvals(A,B); bb,bp=bvals(A,B)
    ok=(cp==bp)
    for k,v in bb.items():
        cv=cb.get(k)
        if cv is None:
            # chomp only prints >=2 and 1x1; re-derive: run chomp with dump? instead trust: if v>=2 must be present
            if v>=2 or k==(1,1): ok=False; print("missing",A,B,k,v)
            continue
        if cv!=v: ok=False; print("mismatch",A,B,k,"C",cv,"brute",v)
    print("box %dx%d: Pcount C=%s brute=%s -> %s"%(A,B,cp,bp,"MATCH" if ok else "FAIL"))
    allok&=ok
print("VALIDATION:","PASS" if allok else "FAIL")

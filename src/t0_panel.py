import sys; sys.path.insert(0,"src")
import numpy as np, requests, proto_tools as pt
import helpers as H

PDB_IDS=["1UBQ","1VII","1ENH","1PGB","1SHG","1CSP","1FKB","1TEN","2GB1","1LMB"]
AA=set("ACDEFGHIKLMNPQRSTVWY"); HYDRO=set("LIVFMACWY")

def fetch_seq(pid):
    t=requests.get(f"https://www.rcsb.org/fasta/entry/{pid}",timeout=30).text
    chains,cur=[],[]
    for ln in t.splitlines():
        if ln.startswith(">"):
            if cur: chains.append("".join(cur)); cur=[]
        else: cur.append(ln.strip())
    if cur: chains.append("".join(cur))
    cand=[c for c in chains if set(c)<=AA and 30<=len(c)<=140]
    return max(cand,key=len) if cand else None

def fold_many(seqs):
    cx=[pt.Complex(chains=[pt.Chain(sequence=s)]) for s in seqs]
    out=pt.run_esmfold(pt.ESMFoldInput(complexes=cx),pt.ESMFoldConfig(device="modal"))
    assert len(out.structures)==len(seqs),(len(out.structures),len(seqs))
    return [(st.metrics.model_dump()["avg_plddt"],st.structure) for st in out.structures]

# --- fetch + validate ---
seqs={}
for pid in PDB_IDS:
    try:
        s=fetch_seq(pid); 
        if s: seqs[pid]=s
        print(f"{pid}: len={len(s) if s else 'DROP (no valid single protein chain)'}")
    except Exception as e: print(f"{pid}: fetch error {e}")
names=list(seqs)
print(f"\nusing {len(names)} proteins: {names}\n")

# --- round 1: fold WTs (one batched call) ---
wt=dict(zip(names,fold_many([seqs[n] for n in names])))

# --- build mutants at hydrophobic buried positions only ---
def make_mutants(name):
    seq=seqs[name]; p0,pdb0=wt[name]
    buried=[p for p in H.buried_positions(pdb0) if seq[p] in HYDRO]
    out={}
    for k in (3,6,9):
        if len(buried)>=k:
            m=list(seq)
            for p in buried[:k]: m[p]="D"
            out[k]="".join(m)
    return buried,out

jobs=[]  # (name,k,seq)
meta={}
for n in names:
    buried,muts=make_mutants(n); meta[n]=(buried,muts)
    for k,ms in muts.items(): jobs.append((n,k,ms))

# --- round 2: fold all mutants (one batched call) ---
mres=dict(zip([(n,k) for n,k,_ in jobs], fold_many([s for _,_,s in jobs])))

# --- report ---
print(f"{'protein':8} {'len':>4} {'#core':>5}  WT_pLDDT   " + "  ".join(f"k={k}:pLDDT/TM/verdict" for k in (3,6,9)))
viol_count={3:0,6:0,9:0}; tot={3:0,6:0,9:0}
for n in names:
    p0,pdb0=wt[n]; buried,muts=meta[n]
    row=f"{n:8} {len(seqs[n]):>4} {len(buried):>5}  {p0:.2f}      "
    for k in (3,6,9):
        if k in muts:
            p1,pdb1=mres[(n,k)]; tm=H.tmscore(pdb0,pdb1)
            v=(tm>=0.5) and (p1>=0.7*p0); tot[k]+=1; viol_count[k]+= 1 if v else 0
            row+=f"  {p1:.2f}/{tm:.2f}/{'VIOL' if v else 'resp'}"
        else: row+="  --/--/n.a."
    print(row)
print("\nViolation rate (fold+confidence retained despite destabilizing core mutations):")
for k in (3,6,9):
    if tot[k]: print(f"  k={k}: {viol_count[k]}/{tot[k]} proteins")

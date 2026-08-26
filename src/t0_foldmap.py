"""T0 fold-map: characterize WHERE ESMFold is invariant to destabilizing core mutations,
using structural descriptors COMPUTED from each predicted structure (no hand-labeled folds).
Invariance score = TM(mutant_k3, WT); higher = model more blind. Batched Modal calls."""
import sys; sys.path.insert(0,"src")
import numpy as np, requests, proto_tools as pt
import biotite.structure as struc
import helpers as H

IDS=["1UBQ","1VII","1ENH","1PGB","1SHG","1CSP","1FKB","1TEN","2GB1","1LMB",
     "1R69","1BDD","256B","2CI2","1CTF","3CHY","1URN","1TIT","1MJC","1SRL",
     "2PTL","1WIT","1AYE","2CRO","1PRB","1HRC","1E0L","2ACY","1OPD","1RIS",
     "1PGX","1CSK","1NYF","1SHF","1FYN","1BF4","1DIV","1POH","1YCC","1PSF"]
AA=set("ACDEFGHIKLMNPQRSTVWY"); HYDRO=set("LIVFMACWY")

def fetch_seq(pid):
    try: t=requests.get(f"https://www.rcsb.org/fasta/entry/{pid}",timeout=25).text
    except Exception: return None
    chains,cur=[],[]
    for ln in t.splitlines():
        if ln.startswith(">"):
            if cur: chains.append("".join(cur)); cur=[]
        else: cur.append(ln.strip())
    if cur: chains.append("".join(cur))
    cand=[c for c in chains if set(c)<=AA and 30<=len(c)<=140]
    return max(cand,key=len) if cand else None

def fold_many(seqs, chunk=20):
    res=[]
    for i in range(0,len(seqs),chunk):
        part=seqs[i:i+chunk]
        cx=[pt.Complex(chains=[pt.Chain(sequence=s)]) for s in part]
        out=pt.run_esmfold(pt.ESMFoldInput(complexes=cx),pt.ESMFoldConfig(device="modal"))
        assert len(out.structures)==len(part)
        res+=[(st.metrics.model_dump()["avg_plddt"],st.structure) for st in out.structures]
    return res

def sse_frac(pdb):
    arr=H._load(pdb)
    try: sse=struc.annotate_sse(arr)
    except Exception: return (np.nan,np.nan)
    if len(sse)==0: return (np.nan,np.nan)
    return (float(np.mean(sse=='a')), float(np.mean(sse=='b')))  # helix, sheet fraction

# fetch + validate
seqs={}
for pid in IDS:
    s=fetch_seq(pid)
    if s: seqs[pid]=s
names=list(seqs); print(f"validated {len(names)}/{len(IDS)} proteins")

# WT fold
wt=dict(zip(names,fold_many([seqs[n] for n in names])))

# build k=3 mutants at buried hydrophobic core
rows=[]; mut_seqs=[]; mut_names=[]
for n in names:
    seq=seqs[n]; p0,pdb0=wt[n]
    buried=[p for p in H.buried_positions(pdb0) if seq[p] in HYDRO]
    hf,sf=sse_frac(pdb0)
    rows.append(dict(name=n,length=len(seq),plddt=p0,helix=hf,sheet=sf,ncore=len(buried)))
    if len(buried)>=3:
        m=list(seq)
        for p in buried[:3]: m[p]="D"
        mut_seqs.append("".join(m)); mut_names.append(n)

# mutant fold
mres=dict(zip(mut_names,fold_many(mut_seqs)))
for r in rows:
    n=r["name"]
    if n in mres:
        p1,pdb1=mres[n]; r["tm_k3"]=H.tmscore(wt[n][1],pdb1); r["plddt_mut"]=p1
        r["viol"]=(r["tm_k3"]>=0.5) and (p1>=0.7*r["plddt"])
    else:
        r["tm_k3"]=np.nan; r["viol"]=None

# report
import statistics as st
tested=[r for r in rows if not np.isnan(r["tm_k3"])]
print(f"\n{'protein':7}{'len':>5}{'helix%':>8}{'sheet%':>8}{'ncore':>6}{'WTpLDDT':>9}{'TM_k3':>7}  verdict")
for r in sorted(tested,key=lambda x:-x["tm_k3"]):
    print(f"{r['name']:7}{r['length']:>5}{r['helix']*100:>7.0f}%{r['sheet']*100:>7.0f}%{r['ncore']:>6}{r['plddt']:>9.2f}{r['tm_k3']:>7.2f}  {'VIOL' if r['viol'] else 'resp'}")

# analysis: does invariance (TM_k3) correlate with structural descriptors?
def corr(a,b):
    a,b=np.array(a),np.array(b); m=~(np.isnan(a)|np.isnan(b))
    return float(np.corrcoef(a[m],b[m])[0,1]) if m.sum()>2 else float('nan')
tm=[r["tm_k3"] for r in tested]
print("\n== Pearson r of invariance (TM_k3) vs descriptor ==")
for key in ("helix","sheet","length","ncore","plddt"):
    print(f"   {key:7}: r = {corr([r[key] for r in tested],tm):+.2f}")

# group by dominant secondary structure (COMPUTED, not hand-labeled)
def cls(r):
    if np.isnan(r['helix']): return "?"
    if r['helix']>=0.45 and r['helix']>r['sheet']: return "alpha"
    if r['sheet']>=0.30 and r['sheet']>=r['helix']: return "beta"
    return "mixed"
print("\n== violation rate by computed structure class ==")
from collections import defaultdict
g=defaultdict(list)
for r in tested: g[cls(r)].append(r)
for c,rs in sorted(g.items()):
    vr=sum(1 for r in rs if r["viol"])/len(rs)
    print(f"   {c:6}: n={len(rs):2}  mean TM_k3={np.mean([r['tm_k3'] for r in rs]):.2f}  violation_rate={vr:.0%}")

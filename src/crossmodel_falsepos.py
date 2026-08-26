"""Characterize the cross-model relation's FALSE-POSITIVE regime: for real (natural) proteins,
when do ESMFold and Boltz-2 DISAGREE (low cross-TM)? If disagreement is common for some class of
real proteins, cross-model disagreement is a weak trust filter there. Correlate cross-TM with
size, secondary-structure content, and confidence.
"""
import sys
sys.path.insert(0, "src")
import numpy as np, requests, proto_tools as pt
import biotite.structure as struc
import helpers as H
from scipy import stats

AA = set("ACDEFGHIKLMNPQRSTVWY")
IDS = ["1VII", "1ENH", "1PGB", "1SHG", "1SHF", "1FYN", "1CSP", "1UBQ", "2CI2",
       "1TEN", "1LMB", "1URN", "1RIS", "256B", "1FKB", "3CHY"]

def fetch(pid):
    t = requests.get(f"https://www.rcsb.org/fasta/entry/{pid}", timeout=25).text
    ch, cur = [], []
    for ln in t.splitlines():
        if ln.startswith(">"):
            if cur: ch.append("".join(cur)); cur = []
        else: cur.append(ln.strip())
    if cur: ch.append("".join(cur))
    c = [s for s in ch if set(s) <= AA and 30 <= len(s) <= 135]
    return max(c, key=len) if c else None

def esm(seq):
    o = pt.run_esmfold(pt.ESMFoldInput(complexes=[pt.Complex(chains=[pt.Chain(sequence=seq)])]),
                       pt.ESMFoldConfig(device="modal"))
    return o.structures[0].metrics.model_dump()["avg_plddt"], o.structures[0].structure

def boltz(seq):
    o = pt.run_boltz2(pt.Boltz2Input(complexes=[pt.Complex(chains=[pt.Chain(sequence=seq)])]),
                      pt.Boltz2Config(device="modal", use_msa=False))
    return o.structures[0].metrics.model_dump().get("confidence_score"), o.structures[0].structure

def sse(pdb):
    arr = H._load(pdb)
    try: s = struc.annotate_sse(arr)
    except Exception: return (np.nan, np.nan)
    return (float(np.mean(s == 'a')), float(np.mean(s == 'b'))) if len(s) else (np.nan, np.nan)

seqs = {p: s for p in IDS if (s := fetch(p))}
rows = []
print(f"{'prot':6}{'len':>5}{'helix':>7}{'sheet':>7}{'ESM':>6}{'Boltz':>7}{'crossTM':>9}", flush=True)
for n, s in seqs.items():
    ep, epdb = esm(s); bc, bpdb = boltz(s)
    h, b = sse(epdb); tm = H.tmscore(epdb, bpdb)
    rows.append(dict(prot=n, length=len(s), helix=h, sheet=b, esm=ep, boltz=bc, cross_tm=tm))
    print(f"{n:6}{len(s):>5}{h*100:>6.0f}%{b*100:>6.0f}%{ep:>6.2f}{bc:>7.2f}{tm:>9.2f}", flush=True)

tm = np.array([r["cross_tm"] for r in rows])
print(f"\nreal-protein cross-model agreement: mean={tm.mean():.2f} min={tm.min():.2f} "
      f"; #disagree(<0.7)={int((tm<0.7).sum())}/{len(tm)}", flush=True)
def corr(key):
    x = np.array([r[key] for r in rows], float); m = ~np.isnan(x)
    r, p = stats.pearsonr(x[m], tm[m]); return r, p
for k in ("length", "helix", "sheet", "esm", "boltz"):
    r, p = corr(k); print(f"  cross_tm vs {k:7}: r={r:+.2f} p={p:.3f}", flush=True)
print("\nREAD: which real proteins disagree, and does disagreement track a class (small? beta? "
      "low-confidence?)? That defines where cross-model agreement is/ISN'T a reliable trust filter.", flush=True)
import json
json.dump(rows, open("results/crossmodel_falsepos_result.json", "w"), indent=2, default=str)

"""CRITICAL control: is ESMFold ignoring DESTABILIZING mutations specifically, or is it just
insensitive to ALL mutations? For each protein make two matched k-mutant variants:
  - buried  : k buried-core hydrophobic residues -> Asp  (catastrophic; SHOULD unfold)
  - surface : k most-exposed residues            -> Asp  (benign; should NOT change fold)
Same substitution (->Asp), same count, only LOCATION differs. Discriminator per protein:
  delta = TM(surface,WT) - TM(buried,WT).
A faithful model responds MORE to buried (delta large positive). If delta ~ 0, the model treats a
catastrophic core change like a benign surface one -> it genuinely cannot detect destabilization
(this makes our 'violation' finding specific, not just general mutation-insensitivity).
"""
import sys
sys.path.insert(0, "src")
import numpy as np, requests, proto_tools as pt
import biotite.structure as struc
import helpers as H
from scipy import stats

AA = set("ACDEFGHIKLMNPQRSTVWY"); HYDRO = set("LIVFMACWY")
IDS = ["1UBQ", "1FKB", "1TEN", "1CSP", "1R69", "256B", "1URN", "3CHY", "2ACY", "1YCC"]

def fetch(pid):
    t = requests.get(f"https://www.rcsb.org/fasta/entry/{pid}", timeout=25).text
    ch, cur = [], []
    for ln in t.splitlines():
        if ln.startswith(">"):
            if cur: ch.append("".join(cur)); cur = []
        else: cur.append(ln.strip())
    if cur: ch.append("".join(cur))
    c = [s for s in ch if set(s) <= AA and 60 <= len(s) <= 130]
    return max(c, key=len) if c else None

def fold_many(seqs, chunk=20):
    res = []
    for i in range(0, len(seqs), chunk):
        part = seqs[i:i + chunk]
        o = pt.run_esmfold(pt.ESMFoldInput(complexes=[pt.Complex(chains=[pt.Chain(sequence=s)]) for s in part]),
                           pt.ESMFoldConfig(device="modal"))
        res += [(st.metrics.model_dump()["avg_plddt"], st.structure) for st in o.structures]
    return res

def rsasa(pdb):
    arr = H._load(pdb); a = struc.sasa(arr, vdw_radii="Single")
    rs = struc.apply_residue_wise(arr, a, np.nansum)
    _, names = struc.get_residues(arr)
    return np.array([rs[i] / H.MAX_ASA.get(H.T2O.get(n, 'A'), 129) for i, n in enumerate(names)])

seqs = {p: s for p in IDS if (s := fetch(p))}
names = list(seqs)
wt = dict(zip(names, fold_many([seqs[n] for n in names])))

K = 3
jobs, meta = [], {}
for n in names:
    s = seqs[n]; p0, pdb0 = wt[n]; rs = rsasa(pdb0)
    buried = [i for i in np.argsort(rs) if s[i] in HYDRO and rs[i] < 0.15][:K]          # most buried hydrophobic
    surface = [int(i) for i in np.argsort(rs)[::-1] if rs[i] > 0.40][:K]                 # most exposed
    if len(buried) < K or len(surface) < K:
        continue
    bmut = list(s); [bmut.__setitem__(i, "D") for i in buried]; bmut = "".join(bmut)
    smut = list(s); [smut.__setitem__(i, "D") for i in surface]; smut = "".join(smut)
    meta[n] = dict(buried=buried, surface=surface)
    jobs += [(n, "buried", bmut), (n, "surface", smut)]

mres = dict(zip([(n, k) for n, k, _ in jobs], fold_many([q for _, _, q in jobs])))
print(f"{'protein':7} {'TM_buried':>9} {'TM_surface':>10} {'delta(surf-bur)':>15}", flush=True)
deltas, tmb_all, tms_all = [], [], []
for n in meta:
    _, pdb0 = wt[n]
    pb, pdbb = mres[(n, "buried")]; ps, pdbs = mres[(n, "surface")]
    tmb = H.tmscore(pdb0, pdbb); tms = H.tmscore(pdb0, pdbs)
    d = tms - tmb; deltas.append(d); tmb_all.append(tmb); tms_all.append(tms)
    print(f"{n:7} {tmb:>9.3f} {tms:>10.3f} {d:>15.3f}", flush=True)

import statistics as st
print(f"\nmean TM_buried (catastrophic) = {st.mean(tmb_all):.3f}", flush=True)
print(f"mean TM_surface (benign)      = {st.mean(tms_all):.3f}", flush=True)
print(f"mean delta (surface - buried) = {st.mean(deltas):+.3f}", flush=True)
if len(deltas) >= 5:
    w, p = stats.wilcoxon(tms_all, tmb_all, alternative="greater")
    print(f"Wilcoxon (TM_surface > TM_buried): p={p:.3f}", flush=True)
    print("READ: delta~0 / large p => ESMFold treats catastrophic core change like a benign surface "
          "change (cannot detect destabilization; violation is SPECIFIC, not general insensitivity). "
          "delta large / small p => ESMFold DOES respond more to destabilization (weakens the claim).", flush=True)
import json
json.dump(dict(names=list(meta), tm_buried=tmb_all, tm_surface=tms_all, deltas=deltas),
          open("results/surface_control_result.json", "w"),
          indent=2)

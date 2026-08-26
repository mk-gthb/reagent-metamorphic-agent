"""Close out the destabilization question: does ESMFold's extra response to buried-destabilizing
(vs benign surface) mutations GROW with the number of mutations k? For larger proteins, sweep
k in {3,6,9,12}, buried-core->Asp vs most-exposed->Asp, measure TM to WT.
  delta(k) = mean TM_surface(k) - mean TM_buried(k)
If delta grows with k, ESMFold DOES eventually detect destabilization (sharpens claim).
If delta stays ~0 even at k=12, it is genuinely insensitive to severe destabilization.
"""
import sys
sys.path.insert(0, "src")
import numpy as np, requests, proto_tools as pt
import biotite.structure as struc
import helpers as H

AA = set("ACDEFGHIKLMNPQRSTVWY"); HYDRO = set("LIVFMACWY")
IDS = ["1FKB", "3CHY", "1URN", "1YCC", "256B", "2ACY", "1TEN", "1RIS"]
KS = [3, 6, 9, 12]

def fetch(pid):
    t = requests.get(f"https://www.rcsb.org/fasta/entry/{pid}", timeout=25).text
    ch, cur = [], []
    for ln in t.splitlines():
        if ln.startswith(">"):
            if cur: ch.append("".join(cur)); cur = []
        else: cur.append(ln.strip())
    if cur: ch.append("".join(cur))
    c = [s for s in ch if set(s) <= AA and 85 <= len(s) <= 135]
    return max(c, key=len) if c else None

def fold_many(seqs, chunk=24):
    r = []
    for i in range(0, len(seqs), chunk):
        part = seqs[i:i + chunk]
        o = pt.run_esmfold(pt.ESMFoldInput(complexes=[pt.Complex(chains=[pt.Chain(sequence=s)]) for s in part]),
                           pt.ESMFoldConfig(device="modal"))
        r += [st.structure for st in o.structures]
    return r

def rsasa(pdb):
    arr = H._load(pdb); a = struc.sasa(arr, vdw_radii="Single")
    rs = struc.apply_residue_wise(arr, a, np.nansum); _, nm = struc.get_residues(arr)
    return np.array([rs[i] / H.MAX_ASA.get(H.T2O.get(n, 'A'), 129) for i, n in enumerate(nm)])

def mut(seq, pos):
    s = list(seq); [s.__setitem__(i, "D") for i in pos]; return "".join(s)

seqs = {p: s for p in IDS if (s := fetch(p))}
names = list(seqs)
wtpdb = dict(zip(names, fold_many([seqs[n] for n in names])))

jobs, meta = [], {}
for n in names:
    s = seqs[n]; rs = rsasa(wtpdb[n])
    buried = [int(i) for i in np.argsort(rs) if s[i] in HYDRO and rs[i] < 0.15]
    surface = [int(i) for i in np.argsort(rs)[::-1] if rs[i] > 0.40]
    if len(buried) < max(KS) or len(surface) < max(KS):
        continue
    meta[n] = (buried, surface)
    for k in KS:
        jobs.append((n, "buried", k, mut(s, buried[:k])))
        jobs.append((n, "surface", k, mut(s, surface[:k])))

folded = dict(zip([(n, t, k) for n, t, k, _ in jobs], fold_many([q for _, _, _, q in jobs])))
print(f"proteins used: {list(meta)}\n", flush=True)
print(f"{'k':>3} {'TM_buried':>10} {'TM_surface':>11} {'delta':>8}", flush=True)
for k in KS:
    tmb = [H.tmscore(wtpdb[n], folded[(n, "buried", k)]) for n in meta]
    tms = [H.tmscore(wtpdb[n], folded[(n, "surface", k)]) for n in meta]
    import statistics as st
    print(f"{k:>3} {st.mean(tmb):>10.3f} {st.mean(tms):>11.3f} {st.mean(tms)-st.mean(tmb):>+8.3f}", flush=True)
print("\nREAD: if delta grows with k, ESMFold detects severe destabilization eventually; "
      "if delta stays ~0 (both TM high), it is insensitive even to gutting the whole core.", flush=True)

"""Strengthen R9/R8 from n=1 to a distribution: optimize N independent sequences purely for
ESMFold pLDDT, then cross-check ALL of them + a set of real proteins against Boltz-2.
Defensible claim if it holds: optimized seqs systematically get high ESMFold confidence that
Boltz-2 does NOT corroborate, while real proteins get cross-model agreement.
Honest scope: this tests uncorroborated confidence, not 'objective nonsense'; and it is asymmetric
(we optimized against ESMFold). A symmetric test (optimize vs Boltz, check ESMFold) is future work.
"""
import sys, random, math, json
sys.path.insert(0, "src")
import requests, proto_tools as pt
import helpers as H

AA = "ACDEFGHIKLMNPQRSTVWY"

def esm_batch(seqs):
    o = pt.run_esmfold(pt.ESMFoldInput(complexes=[pt.Complex(chains=[pt.Chain(sequence=s)]) for s in seqs]),
                       pt.ESMFoldConfig(device="modal"))
    return [(st.metrics.model_dump()["avg_plddt"], st.metrics.model_dump().get("ptm"), st.structure)
            for st in o.structures]

def boltz_one(seq):
    o = pt.run_boltz2(pt.Boltz2Input(complexes=[pt.Complex(chains=[pt.Chain(sequence=seq)])]),
                      pt.Boltz2Config(device="modal", use_msa=False))
    m = o.structures[0].metrics.model_dump()
    return m.get("confidence_score"), m.get("complex_plddt"), m.get("ptm"), o.structures[0].structure

def mutate(seq, rng):
    s = list(seq); s[rng.randrange(len(s))] = rng.choice(AA); return "".join(s)

def hillclimb(seed, L=60, rounds=12, batch=10):
    rng = random.Random(seed)
    cur = "".join(rng.choice(AA) for _ in range(L))
    p, ptm, pdb = esm_batch([cur])[0]
    for r in range(1, rounds + 1):
        props = [mutate(cur, rng) for _ in range(batch)]
        res = esm_batch(props)
        bi = max(range(batch), key=lambda i: res[i][0])
        T = max(0.02, 0.3 * (1 - r / rounds))
        if res[bi][0] > p or rng.random() < math.exp((res[bi][0] - p) / T):
            cur, (p, ptm, pdb) = props[bi], res[bi]
    print(f"  seed {seed}: ESMFold pLDDT={p:.3f} pTM={ptm:.3f}", flush=True)
    return cur, p, ptm, pdb

def fetch(pid):
    t = requests.get(f"https://www.rcsb.org/fasta/entry/{pid}", timeout=25).text
    ch, cur = [], []
    for ln in t.splitlines():
        if ln.startswith(">"):
            if cur: ch.append("".join(cur)); cur = []
        else: cur.append(ln.strip())
    if cur: ch.append("".join(cur))
    S = set(AA); return max((c for c in ch if set(c) <= S and 30 <= len(c) <= 140), key=len)

print("optimizing 5 sequences...", flush=True)
opt = [hillclimb(s) for s in range(5)]
reals = {pid: fetch(pid) for pid in ["1UBQ", "1FKB", "1SHG", "1ENH", "1CSP"]}
real_esm = dict(zip(reals, esm_batch(list(reals.values()))))

def line(tag, esm_plddt, esm_pdb, seq):
    bc, bpl, bptm, bpdb = boltz_one(seq)
    tm = H.tmscore(esm_pdb, bpdb)
    gap = esm_plddt - (bc if bc is not None else float('nan'))
    print(f"{tag:16} ESMpLDDT={esm_plddt:.2f}  BoltzConf={bc:.2f}  gap={gap:+.2f}  crossTM={tm:.2f}", flush=True)
    return dict(tag=tag, esm_plddt=esm_plddt, boltz_conf=bc, gap=gap, cross_tm=tm)

print("\n=== OPTIMIZED (optimized against ESMFold) ===", flush=True)
rows = [line(f"opt_seed{i}", p, pdb, seq) for i, (seq, p, ptm, pdb) in enumerate(opt)]
print("\n=== REAL proteins (controls) ===", flush=True)
rows += [line(f"real_{pid}", real_esm[pid][0], real_esm[pid][2], reals[pid]) for pid in reals]

import statistics as st
og = [r["gap"] for r in rows if r["tag"].startswith("opt")]
rg = [r["gap"] for r in rows if r["tag"].startswith("real")]
oc = [r["cross_tm"] for r in rows if r["tag"].startswith("opt")]
rc = [r["cross_tm"] for r in rows if r["tag"].startswith("real")]
print("\n=== SUMMARY ===", flush=True)
print(f"optimized:  mean ESM-Boltz gap={st.mean(og):+.2f}  mean crossTM={st.mean(oc):.2f}  (n={len(og)})")
print(f"real:       mean ESM-Boltz gap={st.mean(rg):+.2f}  mean crossTM={st.mean(rc):.2f}  (n={len(rg)})")
json.dump(rows, open("results/r9_multi_result.json", "w"),
          indent=2, default=str)
print("READ: if optimized has large positive gap + low crossTM while real has ~0 gap + high crossTM, "
      "the pattern holds (ESMFold confidence inflatable without independent corroboration).")

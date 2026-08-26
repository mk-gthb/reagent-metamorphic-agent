"""Size-matched control for the cross-model relation. The n=5 result was confounded: optimized
seqs are all 60 aa and cross-model TM runs low for short chains. So compare the 60-aa optimized
sequences against REAL proteins of the SAME size (~54-66 aa). If optimized STILL show a larger
ESMFold-Boltz gap / lower agreement than size-matched reals, the effect is real. If they overlap,
the cross-model signal was mostly a small-protein artifact (an important honest negative).
"""
import sys, json
sys.path.insert(0, "src")
import requests, proto_tools as pt
import helpers as H
from scipy import stats

AA = set("ACDEFGHIKLMNPQRSTVWY")
opt = json.load(open("results/r9_multi_result.json"))
opt = [r for r in opt if r["tag"].startswith("opt")]  # already have esm_plddt, boltz_conf, gap, cross_tm

def fetch(pid):
    t = requests.get(f"https://www.rcsb.org/fasta/entry/{pid}", timeout=25).text
    ch, cur = [], []
    for ln in t.splitlines():
        if ln.startswith(">"):
            if cur: ch.append("".join(cur)); cur = []
        else: cur.append(ln.strip())
    if cur: ch.append("".join(cur))
    c = [s for s in ch if set(s) <= AA and 54 <= len(s) <= 66]
    return c[0] if c else None

def esm(seq):
    o = pt.run_esmfold(pt.ESMFoldInput(complexes=[pt.Complex(chains=[pt.Chain(sequence=seq)])]),
                       pt.ESMFoldConfig(device="modal"))
    return o.structures[0].metrics.model_dump()["avg_plddt"], o.structures[0].structure

def boltz(seq):
    o = pt.run_boltz2(pt.Boltz2Input(complexes=[pt.Complex(chains=[pt.Chain(sequence=seq)])]),
                      pt.Boltz2Config(device="modal", use_msa=False))
    m = o.structures[0].metrics.model_dump()
    return m.get("confidence_score"), o.structures[0].structure

CAND = ["1PGB", "1SHG", "1BDD", "1SHF", "1FYN", "1ENH", "1E0L", "1NYF", "1CSK", "2GB1", "1R69"]
reals = {}
for pid in CAND:
    s = fetch(pid)
    if s and len(reals) < 8: reals[pid] = s
print("size-matched reals:", {k: len(v) for k, v in reals.items()}, flush=True)

real_rows = []
for pid, s in reals.items():
    ep, epdb = esm(s); bc, bpdb = boltz(s)
    gap = ep - bc; tm = H.tmscore(epdb, bpdb)
    real_rows.append(dict(tag=pid, len=len(s), esm_plddt=ep, boltz_conf=bc, gap=gap, cross_tm=tm))
    print(f"{pid:6} len={len(s)}  ESM={ep:.2f} Boltz={bc:.2f} gap={gap:+.2f} crossTM={tm:.2f}", flush=True)

og = [r["gap"] for r in opt]; oc = [r["cross_tm"] for r in opt]
rg = [r["gap"] for r in real_rows]; rc = [r["cross_tm"] for r in real_rows]
import statistics as st
print("\n=== SIZE-MATCHED COMPARISON (optimized 60aa vs real ~60aa) ===", flush=True)
print(f"optimized (n={len(og)}): gap {st.mean(og):+.2f}  crossTM {st.mean(oc):.2f}")
print(f"real ~60aa (n={len(rg)}): gap {st.mean(rg):+.2f}  crossTM {st.mean(rc):.2f}")
if len(rg) >= 3:
    ug, pg = stats.mannwhitneyu(og, rg, alternative="greater")
    uc, pc = stats.mannwhitneyu(oc, rc, alternative="less")
    print(f"Mann-Whitney: gap(opt>real) p={pg:.3f} ; crossTM(opt<real) p={pc:.3f}")
    print("READ: small p => optimized genuinely differ from SIZE-MATCHED reals (effect is real, not size). "
          "Large p => the cross-model signal was mostly a short-chain artifact (honest negative).")
json.dump(dict(optimized=opt, real=real_rows), open(
    "results/size_matched_result.json", "w"), indent=2, default=str)

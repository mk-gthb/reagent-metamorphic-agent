"""Cross-model check (relation R8 + the real R9 hack test): do ESMFold and Boltz-2 agree?
Runs both models on the R9-optimized sequence + two natural proteins as references.
If Boltz-2 gives the optimized sequence LOW confidence while agreeing with ESMFold on the
naturals, ESMFold's confidence is uncorroborated -> reward-hack demonstrated. If Boltz-2
corroborates it, the sequence is genuinely foldable (not a hack). Honest either way.
"""
import sys, time, json
sys.path.insert(0, "src")
import requests, proto_tools as pt
import helpers as H

UBQ = "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG"
HACKED = json.load(open("results/r9_result.json"))["hacked_seq"]
def fetch(pid):
    t = requests.get(f"https://www.rcsb.org/fasta/entry/{pid}", timeout=25).text
    chains, cur = [], []
    for ln in t.splitlines():
        if ln.startswith(">"):
            if cur: chains.append("".join(cur)); cur = []
        else: cur.append(ln.strip())
    if cur: chains.append("".join(cur))
    AA = set("ACDEFGHIKLMNPQRSTVWY")
    return max((c for c in chains if set(c) <= AA and 30 <= len(c) <= 140), key=len)

def esmfold(seq):
    o = pt.run_esmfold(pt.ESMFoldInput(complexes=[pt.Complex(chains=[pt.Chain(sequence=seq)])]),
                       pt.ESMFoldConfig(device="modal"))
    return o.structures[0].metrics.model_dump(), o.structures[0].structure

def boltz(seq, retries=20, wait=30):
    for i in range(retries):
        try:
            o = pt.run_boltz2(pt.Boltz2Input(complexes=[pt.Complex(chains=[pt.Chain(sequence=seq)])]),
                              pt.Boltz2Config(device="modal", use_msa=False))
            return o.structures[0].metrics.model_dump(), o.structures[0].structure
        except Exception as e:
            print(f"[boltz retry {i+1}/{retries}] {type(e).__name__}: {str(e)[:100]}", flush=True)
            time.sleep(wait)
    raise RuntimeError("Boltz-2 not available after retries")

seqs = {"R9_optimized": HACKED, "ubiquitin(real)": UBQ, "FKBP12(real)": fetch("1FKB")}
print("seqs:", {k: len(v) for k, v in seqs.items()}, flush=True)

rows = {}
for name, s in seqs.items():
    em, epdb = esmfold(s)
    bm, bpdb = boltz(s)
    if name == list(seqs)[0]:
        print("Boltz2 metric keys:", list(bm.keys()), flush=True)
    cross_tm = H.tmscore(epdb, bpdb)  # do the two models agree on the structure?
    rows[name] = dict(esm_plddt=em.get("avg_plddt"), esm_ptm=em.get("ptm"),
                      boltz=bm, cross_tm=cross_tm)
    print(f"\n{name}: ESMFold pLDDT={em.get('avg_plddt'):.3f} pTM={em.get('ptm'):.3f} | "
          f"Boltz2 metrics={ {k:(round(v,3) if isinstance(v,float) else v) for k,v in bm.items()} } | "
          f"cross-model TM(ESM,Boltz)={cross_tm:.3f}", flush=True)

json.dump({k: {kk: vv for kk, vv in r.items()} for k, r in rows.items()},
          open("results/cross_model_result.json", "w"),
          indent=2, default=str)
print("\nDONE. Read: compare Boltz confidence on R9_optimized vs the two real proteins, and the "
      "cross-model TM. Low Boltz confidence / low cross-TM on R9_optimized (but high on reals) = hack.")

"""R9 — is ESMFold pLDDT hackable? A transparent batched hill-climb (with MCMC-style acceptance)
optimizes a RANDOM sequence purely for pLDDT, via Proto's ESMFold on Modal. Oracle-free controls
expose the hack: (a) random baseline hit-rate, (b) pTM (global coherence) vs pLDDT (local confidence).

Uses our own search (not Proto's gradient optimizer) so every step is inspectable. The model under
test is still ESMFold via Proto; only the search loop is ours.
"""
import sys, random, json, math
sys.path.insert(0, "src")
import proto_tools as pt

AA = "ACDEFGHIKLMNPQRSTVWY"


def fold_batch(seqs):
    o = pt.run_esmfold(
        pt.ESMFoldInput(complexes=[pt.Complex(chains=[pt.Chain(sequence=s)]) for s in seqs]),
        pt.ESMFoldConfig(device="modal"))
    return [(st.metrics.model_dump()["avg_plddt"], st.metrics.model_dump().get("ptm"))
            for st in o.structures]


def mutate(seq, rng, nmut=1):
    s = list(seq)
    for _ in range(nmut):
        s[rng.randrange(len(s))] = rng.choice(AA)
    return "".join(s)


def hillclimb(L=60, rounds=20, batch=12, seed=0):
    rng = random.Random(seed)
    cur = "".join(rng.choice(AA) for _ in range(L))
    (cur_p, cur_ptm), = fold_batch([cur])
    traj = [(0, cur_p, cur_ptm)]
    for r in range(1, rounds + 1):
        props = [mutate(cur, rng, 1) for _ in range(batch)]
        res = fold_batch(props)
        bi = max(range(batch), key=lambda i: res[i][0])
        bp, bptm = res[bi]
        T = max(0.02, 0.3 * (1 - r / rounds))
        if bp > cur_p or rng.random() < math.exp((bp - cur_p) / T):
            cur, cur_p, cur_ptm = props[bi], bp, bptm
        traj.append((r, cur_p, cur_ptm))
        print(f"round {r:2}: best_pLDDT={cur_p:.3f}  ptm={cur_ptm:.3f}", flush=True)
    return cur, cur_p, cur_ptm, traj


def random_baseline(L=60, n=40, seed=1):
    rng = random.Random(seed)
    seqs = ["".join(rng.choice(AA) for _ in range(L)) for _ in range(n)]
    ps = [p for p, _ in fold_batch(seqs)]
    return sum(1 for p in ps if p > 0.7) / n, max(ps), sum(ps) / n


if __name__ == "__main__":
    seq, p, ptm, traj = hillclimb()
    frac, mx, mean = random_baseline()
    json.dump(dict(hacked_seq=seq, hacked_plddt=p, hacked_ptm=ptm, trajectory=traj,
                   baseline_frac_gt0p7=frac, baseline_max=mx, baseline_mean=mean),
              open("results/r9_result.json", "w"), indent=2)
    print("\n=== R9 SUMMARY ===")
    print(f"random start pLDDT = {traj[0][1]:.3f}")
    print(f"hacked pLDDT       = {p:.3f}   (pTM = {ptm:.3f})")
    print(f"random baseline    : {frac:.0%} of random 60-mers reach pLDDT>0.7 "
          f"(max {mx:.2f}, mean {mean:.2f})")
    print("READ: optimizing ONLY pLDDT inflates it far above the random baseline. If pTM stays "
          "low while pLDDT is high, the metric is locally hacked without global structural coherence.")

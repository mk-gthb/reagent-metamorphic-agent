"""T0 — metamorphic relation R1 (directional): a destabilizing buried-core mutation MUST
degrade the predicted structure or its confidence. If the model keeps the same fold at high
confidence, that is a VIOLATION (the model ignored biophysics). Oracle-free.

Reports a survival curve over mutation count k across many proteins — R1 is a SOFT relation,
so the statistics are the evidence, never a single-mutant verdict.

TODO after `smoke_test.py`: fill `predict()` with the verified ESMFold call + pLDDT/PDB fields.
"""
import proto_tools as pt
from helpers import buried_positions, tmscore

# ~10 small monomeric proteins spanning folds (fill sequences from RCSB by PDB id).
WT = {
    "1ubq": "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG",
    # 1vii, 1enh, 1pgb, 1shg, 1csp, 1fkb, 1lmb, 1ten ...  (see SPEC.md 9.3)
}
DESTAB = {'L': 'D', 'I': 'D', 'V': 'D', 'F': 'D', 'M': 'D', 'A': 'K',
          'C': 'D', 'W': 'D', 'Y': 'D'}  # bulky-hydrophobic -> charged (destabilizing in expectation)


def mutant(seq, positions, k):
    s = list(seq)
    for p in positions[:k]:
        s[p] = DESTAB.get(s[p], "D")
    return "".join(s)


def predict(seq):
    """Return (avg_pLDDT [0-1], pdb_text) from ESMFold on Modal. VERIFIED against Proto API."""
    out = pt.run_esmfold(
        pt.ESMFoldInput(complexes=[pt.Complex(chains=[pt.Chain(sequence=seq)])]),
        pt.ESMFoldConfig(device="modal"),
    )
    s = out.structures[0]
    return s.metrics.model_dump()["avg_plddt"], s.structure   # avg_plddt is 0-1; ptm/avg_pae also available


def run():
    for name, wt in WT.items():
        p0, pdb0 = predict(wt)
        buried = buried_positions(pdb0)            # use the model's OWN predicted core (self-consistent)
        for k in (3, 6, 9):
            p1, pdb1 = predict(mutant(wt, buried, k))
            tm = tmscore(pdb0, pdb1)
            viol = (tm >= 0.5) and (p1 >= 0.7 * p0)  # fold AND confidence retained
            print(f"{name} k={k}: pLDDT {p0:.2f}->{p1:.2f}  TM={tm:.3f}  "
                  f"{'VIOLATION (ignored destabilization)' if viol else 'ok (responded)'}")


if __name__ == "__main__":
    run()

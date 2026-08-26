# re:AGENT — Metamorphic-testing agent for protein foundation models

**Build spec, v1.** Author: Mana. Anchor: protein structure prediction. Breadth demo: protein–ligand binding.
Governing rule: oracle-free where possible; every relation carries an explicit epistemic status; never overclaim novelty (value = method transfer + artifact, not "models bluff").

> This is the *original plan*. For what actually held up under adversarial controls — including claims we tested and then **retracted** — see [RESULTS.md](../RESULTS.md) ("Bottom line").

---

## 0. One-sentence thesis
An autonomous, **self-improving** agent that generates biophysical **metamorphic relations**, runs them across AlphaFold3 / ESMFold / Boltz-2, and distills the violations into **named failure classes** shipped as a BenchFlow benchmark. It needs no ground-truth structures because it tests *relations between paired predictions*, not absolute correctness.

**The self-improvement thesis (why the loop is legitimate here and not elsewhere):** self-improving algorithms in biology normally fail because their reward is gameable — docking, pLDDT, and the refolding oracle are all ~0.66–0.73 AUC or directly hackable (Garcia/Rocklin; Korbeld/Fürst), so a self-improving loop just learns to game its scorer. A **metamorphic verdict** ("did the model violate relation R?") is the rare reward that is cheap, oracle-free, and hard to game. So this is plausibly the first self-improving bio loop whose reward is *honest by construction*. Pitch line: "self-improvement in biology keeps failing because the reward is gameable — we give it one that isn't."

## 1. Why this is defensible (state it, then defend it)
- Not novel: "structure models bluff / miscalibrated confidence." Documented (Feldman 2026, Garcia/Rocklin 2025, Korbeld/Fürst 2025, Chakravarty/Porter 2024).
- Novel (as a *system + artifact*, w/ absence-of-evidence caveat): importing CheckList behavioral testing + automated red-teaming + quality-diversity search + slice discovery to protein FMs, unified and shipped as a living benchmark. No instance found in Paperclip's bio corpus of an autonomous, self-improving, multi-model, test-generating + failure-characterizing agent.
- The bar it must clear: **≥1 failure class not already in the manual papers.** If we only reproduce "AF3 invariant to core mutation," we failed. Multi-model disagreement classes and fold-class-stratified classes are the most likely sources of a genuinely new finding.

## 2. The metamorphic relations (the intellectual asset)
CheckList taxonomy: MFT = easy known case must pass; INV = label-preserving change → output must NOT move; DIR = change → output must move a specific way.
Each relation below tagged with **epistemic status**: [HARD] = oracle-free/self-contained; [GT] = uses legit experimental ground truth; [SOFT] = relies on a biophysical prior that can be wrong (must aggregate statistically, never claim per-instance truth).

| ID | Type | Relation | Measure | Violation | Status |
|----|------|----------|---------|-----------|--------|
| R1 | DIR | Destabilizing buried-core substitution (e.g. Leu/Ile/Val buried → Asp/Arg) must degrade structure or confidence | TM-score(pred_mut, pred_wt); ΔpLDDT | high TM AND retained pLDDT after k such mutations | SOFT |
| R2 | DIR | Deleting a secondary-structure element must change the fold | TM-score(pred_del, pred_wt) | fold preserved (TM ≥ 0.5) | SOFT→HARD (deletion is structurally unavoidable) |
| R3 | INV | Conservative surface substitution should be ~invariant | TM-score; ΔpLDDT | large structural swing on a benign change (over-sensitivity) | SOFT |
| R4 | MFT | Known fold-switching proteins must admit ≥2 conformations | compare to both experimental PDB folds | model collapses to one fold | GT |
| R5 | MFT | Well-characterized stable protein → high confidence (positive control) | pLDDT vs known | fails an easy case | GT |
| R6 | DIR/neg | Realistic-but-false sequences (shuffled preserving AA composition; 3rd-order Markov) must NOT get high confidence | FP rate at pLDDT>70 | high-confidence junk (baseline ~1/435, Xu/Salzberg) | HARD |
| R7 | INV | Adding/removing evolutionary info (MSA depth, ordering) shouldn't flip a bad design to high confidence | ΔpLDDT with vs w/o MSA | MSA inflates confidence for a non-folding design | SOFT, model-specific (AF2/ColabFold only) |
| R8 | consistency | Independent models agree on real proteins | pairwise TM across AF3/ESMFold/Boltz | systematic disagreement region | FLAG not violation (one may be right — epistemically weaker; use only for slice discovery, never as "error") |
| R9 | reward-hack | Optimizing a sequence to maximize the designability metric (pLDDT/scRMSD) produces high scores on nonsense targets | can optimizer hit pLDDT>85 for scrambled/random backbones? | yes → metric is hackable (Korbeld/Fürst prediction, live) | HARD (self-contained) |

**Critical notes on soft ground truth (R1/R3):** "destabilizing" is a *prior*, not truth — some buried mutations are tolerated. Mitigations: (a) select positions by low relative solvent accessibility on the predicted structure; (b) apply k mutations and read the *distribution/survival-curve* over many proteins (Feldman's approach), never a single-mutation verdict; (c) report as "the model is statistically invariant to changes that are destabilizing in expectation," not "this specific mutant is unstable and the model is wrong." This keeps R1 honest.

## 3. Architecture
Two nested loops (Darwin–Gödel dual-loop pattern, but with an oracle-free reward):
```
╔══ OUTER LOOP: the tester improves ITSELF (self-improvement) ══════════════╗
║  Claude proposes a modification to the tester's own strategy:             ║
║   new relation / perturbation heuristic / target policy / generation code ║
║        ↓                                                                   ║
║   run INNER LOOP on a FIXED held-out protein set                          ║
║        ↓                                                                   ║
║   accept modification only if validated yield+coverage rises significantly ║
║        ↓  (archive of tester variants; sample parent → mutate → evaluate) ║
╚═══════════════════════════════════════════════════════════════════════════╝
        │ (best current tester)
        ▼
── INNER LOOP: the metamorphic search ─────────────────────────────────────
 Paperclip → propose relations from literature (biophysical rules → MRs)
     ↓
 MAP-Elites archive (quality-diversity) ← generate DIVERSE test cases   [Proto optimizers]
     ↓                                     (ESMFold = fast search; Boltz-2/AF verify elites; Modal)
 Run paired predictions → check relation → log violations w/ severity
     ↓
 Slice discovery (ESM-C embeddings + clustering + Claude captions) → named failure classes
     ↓
 BenchFlow environment (the artifact)  +  live R9 reward-hacking demo
```

### 3.1 MAP-Elites archive
- Behavior-descriptor grid (cells): **fold class** (all-α, all-β, α/β, α+β, membrane, IDR) × **perturbation type** (R1–R7) × **model** (AF3, ESMFold, Boltz-2). Optional 4th axis: size bin.
- Cell content: the single test case with max violation severity found for that cell.
- Fitness (severity), e.g. for R1: `retained_TM × retained_pLDDT_fraction` (high = worse bluff).
- Why QD not plain optimization: avoids mode collapse (HARM's documented failure), and **archive coverage over rounds = the self-improvement metric and the "not too concentrated" guarantee.**

### 3.2 Slice discovery
- Embed violating sequences (ESM-C) and/or structural descriptors; fit mixture model / cluster (DOMINO-style); Claude reads each cluster's members → names the class.
- Honest caveat (from lit): entangled embeddings produce incoherent slices. Validate: silhouette + manual spot-check of ≥3 clusters; report method and coherence. Don't present a cluster as a "class" without a checkable shared property.

### 3.3 Self-improvement loop (the genuine one — Darwin–Gödel-style)
The inner loop (§3.1) is a *search*; on its own it's only weakly self-improving. The outer loop is what makes the *tester itself* improve, DGM-style. Distinct from the DGDM drug-discovery paper (which self-improves against an admittedly imperfect docking surrogate), our outer-loop reward is oracle-free.

**What the agent rewrites (the "self"):**
- the **relation library** — invents/edits MRs (e.g. proposes "delete the longest loop" as a new DIR test);
- **perturbation heuristics** — e.g. learns "target highest native-contact positions" beats "random buried";
- **target-selection policy** — which fold classes/proteins to attack next (active testing);
- its own **generation prompt/code** — the Gödel part: modify the code that proposes tests.

**Reward `R_tester` (oracle-free — the key enabler):**
`R_tester = validated_deduplicated_violations + λ · archive_coverage`, evaluated on a **fixed held-out protein set**. Coverage = # MAP-Elites cells newly filled (diversity). A metamorphic violation is self-evidently a violation → no oracle needed.

**Accept/reject (statistical, anti-noise; mirrors DGDM's PAC safeguard):**
Run parent vs. child tester on the same held-out set (paired). Accept the self-modification only if the gain in `R_tester` is significant (paired test) AND effect size clears a threshold — never accept lucky noise. Keep an archive of tester variants (open-ended evolution); sample a parent, mutate, evaluate, retain if better-or-novel.

**Anti-gaming — meta-reward-hacking (CRITICAL, and ironic):**
Our tester is rewarded to "find violations" — so it can reward-hack its *own* objective by inflating trivial/duplicate/invalid violations. The project studies exactly this failure, so we must not commit it. Guards:
- **Validity filter:** a violation counts toward reward only if (a) the perturbation is in-spec (e.g. mutation truly at a buried position), (b) it comes from a HARD or GT relation (SOFT relations are reported, never rewarded), (c) it is deduplicated against the archive.
- **Diversity term** forces spread across cells — can't farm one exploit.
- **Held-out proteins:** self-modifications must *generalize* (cross-protein transfer), echoing DGM's cross-language transfer; overfitting the training proteins is rejected by the held-out reward.
- **Justification gate:** a newly proposed relation must be biophysically justified (Claude cites the rule, Paperclip-grounded) before it can enter the library.

**Self-improvement metric (what the demo shows):** validated held-out yield + coverage vs. outer-loop iteration — the curve should climb, then plateau. Report the plateau honestly (it's an arms race, not magic).

## 4. Compute / feasibility reality check
- **ESMFold**: fast (~secs/protein, single-seq, GPU) → the workhorse for the MAP-Elites search (needs 100s–1000s of evals). Primary model.
- **Boltz-2**: on Modal, warm weights → run on archive *elites* only (confirmation), not every candidate.
- **AF3**: gated/heavy. Realistic plan: AF2 via ColabFold if a path exists; otherwise AF3 is a stretch goal, and R7 (MSA) is AF2/ColabFold-only. Do NOT promise AF3 in the demo; list as "extend to."
- Pattern (mirrors the deep-research finding): cheap model drives the loop, expensive model certifies elites.
- Everything through Proto's tool interface where possible (model runs + optimizer for R1/R2/R9).

## 5. Demo tiers (build in this order; each is shippable)
- **T0 (must, hours):** harness runs one DIR relation (R1) live on ~10 proteins with ESMFold; show a violation + a non-violation. Proves the loop.
- **T1 (the pitch):** MAP-Elites archive filled across fold × perturbation × ≥2 models; heatmap of violation rates by cell. **Deliver ≥1 new class** (fold-stratified or cross-model). This is the make-or-break.
- **T2 (the wow):** R9 live — optimize a sequence to pLDDT>85 in real time, then reveal it also scores high on a scrambled target ⇒ the designability metric is hackable ⇒ any self-improving loop rewarding it is gameable.
- **T2.5 (the self-improving headline):** run ≥3–5 outer-loop iterations and plot validated held-out violation yield + coverage vs. iteration — the curve climbs. This is the deliverable that earns the word "self-improving." Compute caveat: each outer iteration reruns the inner loop, so cap held-out set size (~20–30 proteins) and inner budget; if compute is tight, fewer iterations with an honest "trend, not asymptote" framing beats faking it. Never seed/script the curve — that would be the exact reward-hacking we're calling out.
- **T3 (artifact):** package as a BenchFlow environment; show Paperclip generating one brand-new relation from a paper.

## 6. Risk register (adversarial, per no-hidden-flaws rule)
1. **"You automated Feldman."** → T1 must yield a new class; lean on multi-model (R8→slices) and fold-stratification. If by Sat night no new class appears, pivot the pitch to the *artifact + reward-hacking demo* (T2/T3), which stand alone.
2. **Soft ground truth (R1/R3).** → statistical aggregation + solvent-accessibility position selection + careful claims (see §2). Prefer HARD relations (R2, R6, R9) for headline numbers.
3. **Compute overrun.** → ESMFold-only is a complete T0/T1; Boltz/AF are additive, not required.
4. **Slice incoherence.** → validate + report; fall back to hand-labeled perturbation-type slices if clustering is noisy (still legitimate, just less automated).
5. **R8 misread as "error."** → strictly a reliability flag; never counted as a model being "wrong."
6. **Circularity in R9.** → the scrambled-target control is self-contained (no external truth needed); this is the cleanest relation — consider leading with it.
7. **Novelty overclaim.** → every slide: "first to *transfer* X to protein FMs," never "first to discover models fail."
8. **Meta-reward-hacking (the ironic one).** The self-improving tester can game its own "find violations" reward with trivial/duplicate/invalid violations. → validity filter (HARD/GT only, in-spec, deduped) + diversity term + held-out generalization + justification gate (see §3.3). If we can't show the loop resists gaming its own objective, we cut the outer loop and ship T0–T2 + artifact rather than claim self-improvement we can't defend.
9. **"Self-improving" is thin.** Weak (search-only) self-improvement ≠ the DGM sense. → the outer loop (§3.3) is the real thing; if it doesn't run in time, say "self-improving by design, outer loop is future work," don't overclaim the inner search as self-improvement.

## 7. Open questions to resolve before/at kickoff
- Which fold-switcher set for R4 (Porter lab's curated set?) — get PDB pairs.
- Does Proto expose ESMFold + Boltz-2 as callable tools + an optimizer we can point at pLDDT? (verify in Proto docs/catalog day 1)
- BenchFlow environment schema — minimal viable env wrapping (relation, model, verdict).
- Merge logistics: does this become the trust layer for Denny's EpiGen / Anshu's diagnostics agent?

## 8. References (read in full via Paperclip unless noted)
- Feldman, Brogi, Skolnick. Adversarial mutations in AF/ESMFold. bioRxiv 2026. doi:10.64898/2026.02.25.708002
- Garcia, Dixit, Rocklin. Zero-shot prediction of design success. bioRxiv 2025. doi:10.1101/2025.07.29.667290
- Korbeld, Viliuga, Fürst. Limitations of the refolding pipeline. bioRxiv 2025. doi:10.64898/2025.12.09.693122
- Xu, Salzberg. Reliability of AI-generated protein structures. bioRxiv 2026. doi:10.64898/2026.06.11.731682
- Chakravarty, Lee, Porter. Alternative folds / AF blind spots. PMC11722503 (2024).
- Ribeiro et al. CheckList (MFT/INV/DIR). ACL 2020. arXiv:2005.04118 [web]
- Zhang et al. ML Testing survey (metamorphic). arXiv:1906.10742 [web]
- HARM automated red-teaming. EMNLP 2024. arXiv:2409.16783 [web]
- DOMINO slice discovery (OpenReview) / Distilling Failures as Directions arXiv:2206.14754 [web]
- MAP-Elites / quality-diversity (Mouret & Clune) [web]

---

## 9. Appendix — starter code (T0 + R9)

**Syntax provenance.** VERIFIED from Proto docs: `Segment(sequence=..., sequence_type="protein")` / `Segment(length=L, ...)`; `structure_plddt_constraint` with `function_config={"structure_tool": "esmfold"|"boltz2"|...}` returning **`1.0 − normalized pLDDT`** (so minimizing energy MAXIMIZES pLDDT); `SemigreedyMutationGenerator(SemigreedyMutationGeneratorConfig(position_weighting="uniform"|"entropy"|"plddt", frozen_positions=[...], exclude_current=True))`; `MCMCOptimizer(constructs=, generators=, constraints=, config=MCMCOptimizerConfig(num_steps=, max_temperature=, ...))`, `.run()`, results at `mcmc.constructs[0].joined_sequences`. **VERIFY at kickoff** (docs 404'd): the `Construct` constructor, exact import paths, and how to read the mean pLDDT value + predicted PDB out of a constraint evaluation.

### T0 — R1: does ESMFold ignore destabilizing core mutations? (oracle-free)
```python
# Direct string mutation (more controllable than a generator for R1); Proto only runs ESMFold.
import numpy as np
from proto_language.core import Segment                      # VERIFY path
from proto_language.constraint import Constraint, structure_plddt_constraint  # VERIFY path

WT = {  # ~10 small monomeric proteins spanning folds (all-a / all-b / a-b)
  "1ubq": "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG",
}
# Buried positions: compute rSASA<0.15 from a reference PDB (DSSP/biotite). Placeholder:
BURIED = {"1ubq": [2, 4, 14, 22, 25, 42, 66]}                 # VERIFY via DSSP
DESTAB = {"L":"D","I":"D","V":"D","F":"D","M":"D","A":"K","C":"D","W":"D","Y":"D"}

def mutant(seq, pos, k):
    s = list(seq)
    for p in pos[:k]:
        s[p] = DESTAB.get(s[p], "D")
    return "".join(s)

def predict(seq):                                            # returns (mean_pLDDT, pdb)
    seg = Segment(sequence=seq, sequence_type="protein")
    c = Constraint(inputs=[seg], function=structure_plddt_constraint,
                   function_config={"structure_tool": "esmfold"})
    energy = c.evaluate()                                    # VERIFY eval API
    return (1 - energy) * 100, c.last_structure              # VERIFY attribute

from tmtools import tm_align   # pip install tmtools  (or shell out to USalign)
def tmscore(pdb_a, pdb_b):     # extract CA coords+seq from each, return tm_norm_chain1
    ...                        # VERIFY helper

for name, wt in WT.items():
    p0, s0 = predict(wt)
    for k in (3, 6, 9):
        p1, s1 = predict(mutant(wt, BURIED[name], k))
        tm = tmscore(s0, s1)
        viol = (tm >= 0.5) and (p1 >= 0.7 * p0)             # fold + confidence retained
        print(f"{name} k={k}: pLDDT {p0:.0f}->{p1:.0f} TM={tm:.2f} "
              f"{'VIOLATION (ignored destabilization)' if viol else 'ok (responded)'}")
# Report the SURVIVAL CURVE over k across all proteins, not single verdicts (R1 is SOFT).
```

### R9 — is pLDDT hackable? Optimize ONLY for pLDDT and expose it (oracle-free)
```python
from proto_language.core import Segment, Construct                                   # VERIFY Construct
from proto_language.generator import SemigreedyMutationGenerator, SemigreedyMutationGeneratorConfig
from proto_language.constraint import Constraint, structure_plddt_constraint         # VERIFY path
from proto_language.optimizer import MCMCOptimizer, MCMCOptimizerConfig              # VERIFY path

L = 80
seg = Segment(length=L, sequence_type="protein")            # random init
construct = Construct(segments=[seg])                        # VERIFY constructor
gen = SemigreedyMutationGenerator(SemigreedyMutationGeneratorConfig(position_weighting="plddt"))
plddt_c = Constraint(inputs=[seg], function=structure_plddt_constraint,
                     function_config={"structure_tool": "esmfold"})   # energy = 1 - pLDDT
cfg = MCMCOptimizerConfig(num_results=1, num_steps=300, max_temperature=0.5,
                          min_temperature=0.001, proposals_per_result=1,
                          temperature_schedule="exponential", tracking_interval=1)
mcmc = MCMCOptimizer(constructs=[construct], generators=[gen], constraints=[plddt_c], config=cfg)
mcmc.run()
hacked = mcmc.constructs[0].joined_sequences[0]
# The gotcha: pLDDT climbs to >85 for a sequence optimized for NOTHING but pLDDT.
# Oracle-free controls that expose the hack:
#  (1) baseline hit-rate: fraction of RANDOM L-mers with pLDDT>70 ~ 1/435 (Xu/Salzberg) — compare.
#  (2) naturalness: ESM pseudo-perplexity of `hacked` vs real proteins (should look unnatural).
#  (3) cross-model: run Boltz-2 on `hacked` — confidence collapse / disagreement = hack exposed.
```

### 9.1 Helpers (standard libs — no Proto internals; `pip install biotite tmtools`)
```python
# helpers.py — buried-position detection (rSASA) + TM-score. Independent of Proto.
import io, numpy as np
import biotite.structure as struc
import biotite.structure.io.pdb as pdb
from tmtools import tm_align

# Tien et al. 2013 theoretical max ASA (Å²) for relative-SASA normalization
MAX_ASA = {'A':129,'R':274,'N':195,'D':193,'C':167,'E':223,'Q':225,'G':104,'H':224,
           'I':197,'L':201,'K':236,'M':224,'F':240,'P':159,'S':155,'T':172,'W':285,'Y':263,'V':174}
T2O = {'ALA':'A','ARG':'R','ASN':'N','ASP':'D','CYS':'C','GLU':'E','GLN':'Q','GLY':'G','HIS':'H',
       'ILE':'I','LEU':'L','LYS':'K','MET':'M','PHE':'F','PRO':'P','SER':'S','THR':'T','TRP':'W','TYR':'Y','VAL':'V'}

def _load(pdb_in):
    """Accept a .pdb path or raw PDB text; return amino-acid AtomArray (model 1)."""
    src = io.StringIO(pdb_in) if "\n" in pdb_in else pdb_in
    arr = pdb.PDBFile.read(src).get_structure(model=1)
    return arr[struc.filter_amino_acids(arr)]

def buried_positions(pdb_in, rsasa_cutoff=0.15):
    """0-based residue indices (in structure order) with relative SASA < cutoff."""
    arr = _load(pdb_in)
    atom_sasa = struc.sasa(arr, vdw_radii="Single")               # per-atom SASA (NaN for excluded)
    res_sasa  = struc.apply_residue_wise(arr, atom_sasa, np.nansum)
    _, res_names = struc.get_residues(arr)
    rsasa = np.array([res_sasa[i] / MAX_ASA.get(T2O.get(n,'A'), 129) for i, n in enumerate(res_names)])
    return np.where(rsasa < rsasa_cutoff)[0].tolist()             # buried core
    # CAVEAT: indices are structure-order; assumes 1 chain, no gaps → aligns to seq string index.
    # If residue numbering has gaps, remap via arr.res_id before using as string offsets.

def _ca(pdb_in):
    arr = _load(pdb_in); ca = arr[arr.atom_name == "CA"]
    _, names = struc.get_residues(ca)
    return ca.coord, "".join(T2O.get(n,'X') for n in names)

def tmscore(pdb_a, pdb_b):
    ca_a, sa = _ca(pdb_a); ca_b, sb = _ca(pdb_b)
    return tm_align(ca_a, ca_b, sa, sb).tm_norm_chain1           # 1.0 = identical fold
```

### 9.2 R9 controls (expose the hack; oracle-free)
```python
# random-baseline hit rate + ESM naturalness. `predict_plddt` = the Proto ESMFold call from §T0.
import numpy as np
AA = "ACDEFGHIKLMNPQRSTVWY"

def random_baseline_hitrate(L=80, n=200, thresh=70.0, rng=np.random.default_rng(0)):
    """Fraction of RANDOM L-mers with pLDDT>thresh (expect ~1/435 ≈ 0.2%, Xu/Salzberg)."""
    hits = sum(predict_plddt("".join(rng.choice(list(AA), L)))[0] > thresh for _ in range(n))
    return hits / n

def esm_pseudo_perplexity(seq, model_name="esm2_t33_650M_UR50D"):
    """Higher = less natural. Compare `hacked` vs real proteins. pip install fair-esm."""
    import torch, esm
    model, alphabet = esm.pretrained.__dict__[model_name]()
    bc = alphabet.get_batch_converter(); model.eval()
    _, _, toks = bc([("s", seq)])
    nll = 0.0
    with torch.no_grad():
        for i in range(1, toks.size(1)-1):                       # mask each position
            m = toks.clone(); m[0, i] = alphabet.mask_idx
            lp = torch.log_softmax(model(m)["logits"][0, i], -1)
            nll -= lp[toks[0, i]].item()
    return float(np.exp(nll / len(seq)))
    # NOTE: Proto may expose an ESM-perplexity/naturalness constraint — prefer it if so (one less dep).
```

### 9.3 T0 protein set (proposed — small monomeric, fold-diverse, all have experimental PDBs)
| PDB | Protein | Fold class | ~len |
|-----|---------|-----------|------|
| 1UBQ | ubiquitin | α/β | 76 |
| 1VII | villin headpiece | all-α | 36 |
| 1ENH | engrailed homeodomain | all-α | 54 |
| 1PGB | protein G B1 | α+β | 56 |
| 1SHG | α-spectrin SH3 | all-β | 62 |
| 1CSP | cold-shock protein B | all-β (OB) | 67 |
| 2GB1 | protein G (NMR) | α+β | 56 |
| 1FKB | FKBP12 | α/β | 107 |
| 1LMB | λ-repressor | all-α | 87 |
| 1TEN | tenascin fibronectin-III | all-β | 90 |

Buried positions: run `buried_positions()` on each experimental PDB (fetch from RCSB) **or** on ESMFold's own WT prediction (self-consistent: we destabilize the core of the structure the model itself predicts). Prefer the latter to avoid PDB-numbering remaps.



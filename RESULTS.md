# Results log

## Bottom line (honest synthesis — read this first)

What survived adversarial controls, as of the latest run:

1. **The harness works and is self-correcting.** It repeatedly caught and walked back our OWN overclaims (see the reversals below). That rigor — not any single "flaw" — is the core demonstration of the trust-layer idea.
2. **Destabilization relation (R1): ESMFold is largely FAITHFUL — initial "flaw" retracted.** It correctly ignores benign surface mutations (TM~0.99 at all k) and responds progressively to buried-core destabilization (Δ vs surface grows 0.04→0.30 as k=3→12). It is only weakly sensitive at low k. So "ESMFold ignores destabilization" was an overclaim; metamorphic-invariance is a WEAK red-team signal here.
3. **Cross-model disagreement (R8): the surviving signal.** Sequences optimized purely for ESMFold confidence are NOT corroborated by an independent model (Boltz-2): size-matched vs real ~60aa proteins, gap p=0.015 and agreement p=0.003. Meaning: ESMFold's confidence can be inflated to a point an independent model won't share — so single-model confidence is not self-sufficient. Caveat: ~1/7 real proteins are false positives (characterizing the false-positive regime now).
4. **R9 optimizability:** ESMFold pLDDT is optimizable far above the random baseline, but within-ESMFold pTM rises with it — so a single model can't reveal its own inflation; the independent cross-check is required.
5. **Fold-map:** low-k invariance tracks protein size/confidence, NOT fold type (helix/sheet r≈0). But since low-k invariance is now known to be mostly benign, this is a minor result.

Net: the defensible contribution is (a) a rigorous, self-auditing red-team harness that adjudicates model-trust claims and catches overclaims (including ours), and (b) evidence that ESMFold confidence is inflatable-without-independent-corroboration. The dramatic "ESMFold ignores biophysics" framing did NOT survive scrutiny — reported straight.

---

## T0 — destabilization invariance (metamorphic relation R1), ESMFold via Modal

First end-to-end run. Target: ubiquitin (1UBQ, 76 aa). Buried core positions from ESMFold's own WT prediction (rSASA < 0.15). Destabilizing substitutions: bulky-hydrophobic core → Asp/Lys. Verdict = VIOLATION if fold is retained (TM ≥ 0.5) AND confidence retained (pLDDT ≥ 0.7× WT).

| perturbation | avg_pLDDT | TM(mut, WT) | verdict |
|---|---|---|---|
| WT | 0.90 | — | well-folded (correct) |
| 3 core → Asp | 0.80 | 0.957 | **VIOLATION** — fold unchanged |
| 5 core → Asp | 0.70 | 0.843 | **VIOLATION** — fold unchanged |
| 8 core → Asp | 0.50 | 0.736 | responded (fold/confidence degrade) |

Reading: gutting ubiquitin's hydrophobic core with 3–5 charged residues — which must destabilize/unfold it — leaves ESMFold's predicted structure 84–96% identical at high confidence. The model only begins to react at 8 core mutations. This is the metamorphic-invariance violation (cf. Feldman/Skolnick 2026), reproduced on our own harness.

Caveats (honest): single protein, SOFT relation (destabilization is a biophysical prior, not measured ground truth) — the real evidence is the survival curve across many proteins and fold classes. `avg_plddt` is on a 0–1 scale.

Next: run the 10-protein fold-diverse panel (SPEC §9.3); add R9 (reward-hacking) and cross-model (Boltz-2) relations.

## T0 panel — 10 fold-diverse proteins (batched, 2 Modal calls)

Mutations restricted to buried HYDROPHOBIC core (I/V/L/F/M → Asp) so the destabilization prior is strong. Verdict = VIOLATION if TM(mut,WT) ≥ 0.5 AND pLDDT ≥ 0.7× WT.

| protein | len | WT pLDDT | k=3 (TM) | k=6 | k=9 |
|---|---|---|---|---|---|
| 1UBQ | 76 | 0.86 | VIOL 0.96 | VIOL 0.81 | VIOL 0.93 |
| 1VII | 36 | 0.87 | resp 0.58 | — | — |
| 1ENH | 54 | 0.85 | VIOL 0.88 | resp | — |
| 1PGB | 56 | 0.85 | resp 0.68 | resp | resp |
| 1SHG | 62 | 0.86 | VIOL 0.73 | resp | resp |
| 1CSP | 67 | 0.86 | VIOL 0.92 | resp | resp |
| 1FKB | 107 | 0.90 | VIOL 1.00 | VIOL 0.99 | VIOL 0.95 |
| 1TEN | 90 | 0.89 | VIOL 0.92 | VIOL | resp |
| 2GB1 | 56 | 0.85 | resp 0.68 | resp | resp |
| 1LMB | 92 | 0.88 | VIOL 0.90 | VIOL | resp |

**Violation rate: k=3 → 7/10, k=6 → 4/9, k=9 → 2/8.**

Key observation: insensitivity is strongly protein/fold-specific. FKBP12 (1FKB) retains TM 0.95–1.00 even at 9 core→Asp mutations; protein G (1PGB/2GB1) responds even at k=3. This heterogeneity is the basis for failure-class discovery.

Internal checks passed: tmscore self=1.0, neg-control TM(ubq,peptide)=0.17, 1PGB≡2GB1 identical (determinism), buried positions rSASA<0.12 & hydrophobic. avg_plddt on 0–1 scale.

## T0 fold-map — 37 proteins, invariance vs COMPUTED structural descriptors

k=3 buried-hydrophobic→Asp across 37 validated small monomers (descriptors computed from each ESMFold structure via biotite; no hand-labeled folds). Invariance score = TM(mutant_k3, WT).

- **Overall violation rate at k=3: 30/36 ≈ 83%.**
- **Secondary-structure class does NOT explain invariance:** Pearson r(TM_k3, helix%) = −0.05, r(sheet%) = +0.05. Violation rate 80% (α), 85% (β), 85% (mixed) — indistinguishable. (This *refuted* the small-panel hypothesis that blindness was fold-type-specific.)
- **Invariance correlates with size / packing / confidence:** r(pLDDT)=+0.57, r(core size)=+0.53, r(length)=+0.45. Responders are mostly small and/or low-confidence (1PGB, 2GB1, 1VII, 2PTL, 1PGX, 1BDD).
- Reading: ESMFold is *more* blind to core destabilization for larger, well-packed, high-confidence predictions.

Stats (computed, not asserted; scipy.stats on the observed values): fold-type is NOT a predictor — helix r=−0.05 (p=0.75), sheet r=+0.06 (p=0.74). Size/packing/confidence ARE: pLDDT r=+0.57 (p=0.0004), core size r=+0.53 (p=0.001), length r=+0.45 (p=0.005). Partial correlations show pLDDT and length each predict invariance INDEPENDENTLY: TM vs pLDDT | length r=+0.58 (p<0.001); TM vs length | pLDDT r=+0.46 (p=0.006) — so they are not merely confounded (earlier caveat was too pessimistic).
Caveats (honest): n is really 35 unique proteins — 1PGB and 2GB1 are byte-identical (56 aa), so they are a DUPLICATE in this correlation, not two independent points (they also re-confirm determinism: identical input → identical output). The length correlation may be partly a TM-score short-chain artifact; pLDDT-retention is the cleanest (TM-independent) signal. SOFT destabilization prior throughout.

Next: partial-correlation follow-up; R9 reward-hacking; cross-model (Boltz-2) disagreement.

## R9 — is ESMFold pLDDT hackable? (transparent hill-climb, 20 rounds, L=60)

Optimized a random 60-mer PURELY for ESMFold pLDDT (batched hill-climb, MCMC acceptance).

- pLDDT climbed **0.36 (random) → 0.77**. Random baseline: **0/40 random 60-mers reach pLDDT>0.7** (max 0.55, mean 0.40). So pLDDT is optimizable far above chance.
- **HONEST CATCH (not yet a hack):** pTM rose *in lockstep* (0.25 → 0.68). ESMFold considers the optimized sequence genuinely coherent by both local (pLDDT) and global (pTM) measures. This demonstrates pLDDT is *optimizable* (expected), NOT that it is *fooled*. A reward-hack requires high confidence on a structure that is actually nonsense — not shown here.
- **The real hack test = independent oracle disagreement.** Next: run this optimized ESMFold-high-confidence sequence through Boltz-2. If Boltz-2 gives low confidence / a different fold, ESMFold's confidence is uncorroborated → that is the reward-hack / cross-model relation (R8).

Optimized seq: `PQCKTNPLSNWHTFLFEYKVYFSDHMSVEASMYHVYRTKCRADPKYSMIMDRWMMFVRDC` (pLDDT 0.74, pTM 0.61).

## Cross-model (R8) — ESMFold's confidence is UNCORROBORATED by Boltz-2 (reliability flag)

Ran ESMFold AND independent Boltz-2 (single-sequence) on the R9-optimized sequence + two real proteins.

| sequence | ESMFold pLDDT/pTM | Boltz-2 conf/pLDDT | cross-model TM(ESM,Boltz) |
|---|---|---|---|
| R9 optimized | 0.74 / 0.61 | **0.51 / 0.53** | **0.49** |
| ubiquitin (real) | 0.86 / 0.83 | 0.93 / 0.93 | 0.98 |
| FKBP12 (real) | 0.90 / 0.89 | 0.90 / 0.91 | 0.96 |

- Real proteins: the two independent models AGREE (both confident, cross-TM 0.96–0.98) → trustworthy.
- R9-optimized: ESMFold confident, but Boltz-2 gives LOW confidence (0.51) and a DIFFERENT structure (cross-TM 0.49). ESMFold's confidence is **uncorroborated** by an independent model — a *reliability flag*. This is CONSISTENT WITH a reward-hack but does NOT by itself prove ESMFold (rather than Boltz) is the hallucinating one: cross-model disagreement identifies unreliability, not culpability.
- Key methodological point: single-model confidence (even ESMFold's global pTM) did NOT reveal the unreliability; the independent cross-check did.
- CORRECTION: an earlier version of this note said "reward-hack confirmed" — that was an overclaim (contradicting the SPEC's own "R8 = flag, not violation"). Corrected above.

Honest caveats: n=1 optimized sequence (need a distribution of optimized seqs vs reals to be statistically solid); Boltz run single-sequence (use_msa=False) for a fair comparison to ESMFold and because the optimized seq has no homologs; cross-TM low could in principle reflect a hard/short target, but Boltz's independently-low confidence corroborates it. Boltz metric = confidence_score / complex_plddt.

## Cross-model strengthened to n=5 — pattern holds at group level BUT confounded (honest)

Optimized 5 independent sequences (seeds 0-4) purely for ESMFold pLDDT, cross-checked all vs Boltz-2, plus 5 real proteins as controls. gap = ESMFold pLDDT − Boltz confidence; crossTM = ESMFold-vs-Boltz structural agreement.

| group | mean gap | mean crossTM |
|---|---|---|
| optimized (n=5, all 60 aa) | +0.09 | 0.40 |
| real (n=5) | −0.02 | 0.86 |

**Two problems the bigger sample exposed (do NOT overclaim):**
1. A REAL protein overlaps the optimized group: 1SHG (SH3, 62 aa) has gap +0.11, crossTM 0.49 — indistinguishable from optimized. So "positive gap / low agreement" is NOT a clean hack signature.
2. CONFOUND with size: optimized seqs are all 60 aa; cross-model TM runs lower for short chains regardless. The low agreement may be partly a small-protein artifact, not purely ESMFold bluffing. The confidence gap (+0.09) is also modest.

**Verdict:** the cross-model effect is real at the group level but modest and confounded by protein length; it does NOT cleanly separate optimized from small real proteins. Clean "hack detector" claim does not survive the larger sample.

**Fix (next):** size-matched controls — compare 60-aa optimized sequences against REAL ~60-aa proteins so length cannot explain the gap. Only then is the cross-model relation defensible.

## Size-matched control — the cross-model effect is REAL (not a short-chain artifact)

Compared the 5 optimized 60-aa sequences against 7 REAL proteins of matched size (54-62 aa), both models, single-sequence.

| group | mean gap (ESM−Boltz) | mean cross-model agreement (TM) |
|---|---|---|
| optimized (n=5, 60 aa) | +0.09 | 0.40 |
| real ~60 aa (n=7) | −0.04 | 0.89 |

Mann-Whitney: gap (opt>real) **p=0.015**; agreement (opt<real) **p=0.003**.

**Verdict:** even controlling for size, optimized sequences get significantly higher ESMFold-vs-Boltz confidence gaps and lower cross-model agreement than real proteins of the same length. The earlier n=5 "size confound" worry is RESOLVED — the effect is real. (Being critical strengthened the claim here, vs weakening it at n=5.)

Honest caveat: NOT a clean classifier — 1SHG (real SH3, 62 aa) is a false positive (agreement 0.50, gap +0.12); ~1/7 real proteins also show model disagreement (some real β-folds genuinely split the two models). The cross-model relation flags unreliability with a nonzero false-positive rate; it is a probabilistic signal, not a binary hack-detector.

## Surface-vs-buried control — TEMPERS the destabilization claim (critical)

For 10 proteins, matched k=3 →Asp mutations at buried-core (catastrophic) vs most-exposed (benign) positions; only LOCATION differs.

- mean TM(buried→WT) = **0.956** (should be low if model detected the destabilization)
- mean TM(surface→WT) = **0.988**
- Δ (surface − buried) = **+0.032**, Wilcoxon p=0.014

**Interpretation (revises earlier framing):**
1. BOTH stay near 1.0 → much of the "invariance to destabilization" is really GENERAL insensitivity to a few point mutations, NOT destabilization-specific blindness. Earlier "ESMFold ignores destabilizing mutations" OVERCLAIMED.
2. But a small, significant signal exists: the model responds slightly more to buried-destabilizing than benign-surface mutations (Δ=0.032, p=0.014). Heterogeneous per protein (3CHY Δ=0.14; 1R69/2ACY ≈0).

**Revised honest claim:** a fold-destroying core mutation perturbs ESMFold's prediction almost as little as a benign surface swap (Δ≈0.03 at k=3) — it cannot be relied on to register even severe destabilizing mutations, but this is largely general point-mutation insensitivity, not a clean destabilization-specific failure. This retro-tempers the fold-map/panel "violation" counts (TM≥0.5 is a weak bar that ~any 3 mutations clear).

Limitation: control was at k=3 only; the buried-vs-surface gap may grow at higher k (next: k-sweep). RETRO-CORRECTION applied to the T0/fold-map narrative above.

## k-sweep — destabilization relation LARGELY EXONERATES ESMFold (major reversal)

Buried-core vs surface →Asp, k=3..12, on 8 larger proteins. TM to WT:

| k | TM_buried (destabilizing) | TM_surface (benign) | Δ |
|---|---|---|---|
| 3 | 0.952 | 0.991 | +0.038 |
| 6 | 0.937 | 0.990 | +0.053 |
| 9 | 0.838 | 0.986 | +0.148 |
| 12 | 0.684 | 0.982 | +0.297 |

- Surface (benign) stays ~0.99 at all k → ESMFold correctly ignores benign mutations (NOT generally insensitive).
- Buried (destabilizing) degrades progressively; Δ grows monotonically 0.04→0.30 → ESMFold DOES respond specifically to core destabilization, scaling with severity.

**Reversal:** the original "ESMFold ignores destabilization / is confidently wrong" thesis is NOT supported. ESMFold behaves largely faithfully here — benign→no change (correct), destabilizing→progressive change (correct), only under-sensitive at low k (3-6). The destabilization/metamorphic-invariance relation is therefore a WEAK red-team signal, not a strong flaw.

**Meta-point (the trust layer working):** our own controls (surface-vs-buried, k-sweep) killed our initial "flaw." The harness catches overclaims — including ours. The signal that survived all controls is CROSS-MODEL disagreement (size-matched significant, p=0.003), i.e. optimizing against ESMFold yields sequences an independent model won't corroborate. That — not destabilization-invariance — is the project's real, defensible finding.

## Cross-model false-positive regime (16 real fold/size-diverse proteins)

ESMFold vs Boltz-2 agreement on natural proteins: mean cross-TM **0.90**, **2/16 disagree** (<0.7): 1SHG (0.50), 2CI2 (0.68).
- Disagreement does NOT track size (r=+0.22, p=0.42) or fold type (sheet r=−0.07, n.s.) — retro-confirms the size-matched result was not a size artifact.
- Disagreement DOES track independent-model confidence: cross-TM vs Boltz confidence r=+0.84 (p<0.001). Models disagree mainly where Boltz is itself uncertain.
- Honest limitation: the 2 real false-positives overlap the optimized sequences on gap/agreement, so cross-model disagreement is a RELIABILITY FLAG (≥1 model unsure / they disagree), NOT a clean adversarial-vs-real classifier.

Usable statement: on natural proteins, cross-model agreement is reliable ~14/16 of the time; disagreements coincide with low independent-model confidence, not with size or fold.

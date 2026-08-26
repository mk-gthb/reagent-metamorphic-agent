# re:AGENT demo — a self-auditing red-team harness for protein models

**One line:** We built a harness that stress-tests protein-structure models (ESMFold, Boltz-2) and adjudicates when their confidence can be trusted — and the headline is that **its own controls caught and corrected our overclaims.**

## The story (honest arc)
1. Built a metamorphic-testing harness on real infra (Proto + ESMFold/Boltz-2 on Modal), self-audited (TM self-test=1.0, negative control=0.17, determinism confirmed).
2. First result looked like a flaw: ESMFold seemed to "ignore" destabilizing core mutations.
3. **Our own adversarial controls walked it back.** Surface-vs-buried + a k-sweep showed ESMFold correctly ignores benign surface mutations (TM~0.99) and responds progressively to core destabilization (gap grows Δ 0.04→0.30 as k=3→12). The "flaw" was mostly general low-k insensitivity — ESMFold is largely faithful. (Figure A.)
4. The signal that **survived** every control: sequences optimized purely for ESMFold confidence are not corroborated by an independent model (Boltz-2). Size-matched vs real ~60aa proteins: agreement p=0.003. So single-model confidence is inflatable-without-corroboration. (Figure B.)
5. Honest limit: it's a reliability *flag*, not a clean adversarial detector — 2/16 real proteins are false positives; disagreement tracks the independent model's own confidence (r=0.84), not size or fold.

## Why it matters
The hackathon theme is "results worth trusting." We demonstrate the trust layer by turning it on ourselves: a red-team rigorous enough to **kill its own flashy claim**. The deliverable is the harness + the honest finding, not a manufactured "ESMFold is broken."

## What to show live
- `figures/honest_results.png` (A: k-sweep; B: cross-model).
- `RESULTS.md` "Bottom line" — the reversal trail, every claim with a control or caveat.
- A live cross-model call: optimized seq → ESMFold confident, Boltz not.

## Stack
Claude Code (orchestration) · Proto (design/optimizer + model runner) · ESMFold + Boltz-2 on Modal (H100) · biotite/tmtools (analysis). Paperclip used for the literature grounding behind the metamorphic relations.

## Honest limitations
Destabilization relation is a weak red-team signal (ESMFold mostly faithful). Cross-model is a probabilistic flag with a nonzero false-positive rate. "Destabilizing" is a biophysical prior, not measured ΔΔG. Self-improvement outer loop designed but not built (empirically lower-value once the strong flaw dissolved).

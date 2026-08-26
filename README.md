# reagent-metamorphic-agent

**A red-team evaluation framework for protein structure-prediction models** — probing *when a model's confidence deserves trust*, without any experimental ground truth.

> Built at **re:AGENT** — a hackathon by **Anthropic**, the **Arc Institute**, and **GXL**. San Francisco, 2026.

[![License: MIT](https://img.shields.io/badge/License-MIT-1a1a1a.svg)](LICENSE)

---

Structure predictors such as ESMFold and Boltz-2 increasingly decide which candidates justify weeks of wet-lab work — on the strength of a confidence score they assign themselves. Checking that score normally requires the very experiments it was meant to save.

This framework evaluates that trust with **no ground truth**, using two oracle-free ideas:

- **Metamorphic testing** — perturb an input in a biophysically-constrained way (destroy a protein's buried hydrophobic core, which *must* destabilize it) and verify the prediction moves as it should. Adapted from software testing.
- **Cross-model corroboration** — read two independently-trained models as independent witnesses; where they diverge, confidence is suspect.

Every claim is paired with an adversarial control — the framework is built to falsify its own hypotheses.

**→ Findings and the full experiment log live in [RESULTS.md](RESULTS.md).**

## Repository

```
├── RESULTS.md      experiment log + honest bottom line
├── src/            experiments — metamorphic tests, controls, cross-model checks
├── results/        machine-readable outputs (JSON)
├── figures/
├── docs/           SPEC.md (original plan) · DEMO.md (pitch)
└── requirements.txt
```

## Run

```bash
python3 -m venv .venv && source .venv/bin/activate        # Python 3.10+
pip install -r requirements.txt
pip install git+https://github.com/evo-design/proto-language.git
modal setup && proto-tools deploy --apps esmfold,boltz2 --env proto-env
python smoke_test.py                                       # end-to-end check
```

Models are served via [Proto](https://proto.evodesign.org) on [Modal](https://modal.com); structures compared by TM-score, solvent accessibility via `biotite`, statistics via `scipy`.

## References

- Feldman, Brogi, Skolnick. *Adversarial sequence mutations in AlphaFold and ESMFold.* bioRxiv (2026).
- Garcia, Dixit, Rocklin. *Evaluating zero-shot prediction of protein design success.* bioRxiv (2025).
- Korbeld, Viliuga, Fürst. *Limitations of the refolding pipeline for de novo protein design.* bioRxiv (2025).
- Ribeiro et al. *Beyond Accuracy: Behavioral Testing of NLP Models with CheckList.* ACL (2020).

## Acknowledgements

Built at **re:AGENT** (Anthropic · Arc Institute · GXL). Tooling: Proto (Evo Design), Modal, and Paperclip (GXL) for literature grounding.

## License

[MIT](LICENSE) © 2026 Manasvini Kothari

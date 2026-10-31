# Building a better cancer therapeutic
## Affinity maturation of the trastuzumab scFv against HER2

---

## 1. The problem

HER2 is overproduced in 20–30% of breast cancers and drives tumour growth. Trastuzumab
(Herceptin) is an antibody that grips HER2's extracellular domain and flags the cell for the
immune system. We were given it as a single-chain fragment (scFv) and asked to make it bind HER2
**more tightly**.

- **Parent:** 242-aa scFv, VH – (GGGGS)3 – VL (PDB 1N8Z chains B + A).
- **Number to beat:** 0.82 nM K_D (measured by BLI on Adaptyv's platform).
- **Rules:** ≥5 mutations from the parent; one chain in order VH–(GGGGS)3–VL; ≤260 aa.

## 2. The key decision: maturation, not de-novo design

We already have a binder that works — so this is **affinity maturation**, not inventing a new
protein. That single framing drove every choice:

- We **edit the paratope** (the six CDR loops that actually touch HER2) and **leave the scaffold
  alone**. Framework and the linker hold the shape; changing them mostly breaks the fold without
  improving binding.
- We therefore **excluded the de-novo tools** (RFdiffusion, BoltzGen, BindCraft) entirely — wrong
  job, and too slow for one day.
- We **protected the load-bearing residues** the literature and structure flag as essential
  (W99, R50, R59, Y105 in VH; H226 in VL) and the four conserved cysteines that hold the domains
  together — none were mutated.

## 3. The target, trimmed intelligently

Trastuzumab binds HER2 **subdomain IV** (epitope ~res 557–603). Rather than co-fold against the
full 630-aa ectodomain (slow, noisier), we trimmed HER2 to **domain IV (UniProt P04626 ~510–652,
143 aa)** — faster and a cleaner interface signal.

## 4. Design strategy — a two-track pool of 15

**Track A — literature-grounded (V01–V05).** Built from published, wet-lab-measured wins:
- **AH06** (Moon 2016): CDR-H3 `D102W + M107F + D108A + Y109L` → ~7.4× tighter; the D102W Trp makes
  a new contact with HER2 Ile613.
- **D31E** (Zhou 2025): CDR-H1 edit, 1.5 → 0.2 nM in a related anti-HER2 antibody.
- **S85N + R87T** (Ng 2025): framework edits that retained binding and improved stability.

**Track B — model-proposed (V06–V15).** Each keeps the validated **D102W + D31E anchor** plus 3
model-chosen CDR substitutions, spread across loops for diversity.

## 5. Models used — and why each

| Model | Type | What it contributed |
|---|---|---|
| **ESM2** | general protein language model | zero-shot scan ranking favourable CDR substitutions (the explore step) |
| **LigandMPNN / ProteinMPNN** | structure/inverse-folding | interface-aware "which residues physically fit" |
| **AbLang2** | antibody-specific language model | cross-check that edits are "antibody-native" |
| **Boltz-2** | structure prediction | co-folds each variant + HER2 domain IV (the scoring step) |
| **ipSAE** | interface scoring | ranks the designs by interface confidence |
| **protein-qc** | developability | liability scan (deamidation, cysteines, aggregation) |

The point of three *design* models (general + structural + antibody-specific) is **consensus**: a
mutation that ESM2 and AbLang2 both favour is a far safer bet than one only a single model liked.
AbLang2 agreed on 14 of 30 model picks and pushed back on 10 where a general model had over-reached
— exactly the antibody-specific correction we wanted.

## 6. Rigour / guardrails (verified in code, twice)

Every variant: ≥5 mutations, ≤260 aa, correct scFv format, linker intact, conserved cysteines and
protected contacts untouched, and **no new deamidation (N-[GST]) or isomerization (D-[GST]) motifs**
— our anchor edits actually *remove* two of the parent's liability motifs. A protein-qc liability
scan confirmed no variant is worse than the parent.

## 7. Validation and ranking

We co-folded all 16 (15 variants + parent) against HER2 domain IV with Boltz-2 and ranked by
**ipSAE_min** — the interface-confidence metric the challenge is judged on. The parent is the
reference line (co-folded at ipTM 0.91 / pTM 0.93 / pLDDT 0.91 — a high-confidence, well-docked
complex, as expected for a real antibody).

> See results.md for the full ranked table.

**An honest word on the metric:** ipSAE_min is an interface-confidence *rank*, not an affinity. No
in-silico tool outputs a trustworthy nM K_D for a protein–protein interface — that is a wet-lab
measurement, which is why Adaptyv makes and measures every design and returns K_D ~3 weeks later.
So we use ipSAE to **discard broken designs and pick a diverse set**, not to claim we have "beaten
0.82 nM" on a computer.

## 8. Why this approach is strong

- **Anchored on real data**, then explored with models — not a blind model sweep.
- **Cross-model consensus** (general + structural + antibody-specific) for every model pick.
- **Developability-aware**, not just affinity-chasing.
- **Fully reproducible and open** — every script, input, score and the model patches are on GitHub.
- **Honest about what the numbers mean** — ipSAE to rank, wet lab to measure.

## 9. Limits and next round

One lab round means we spread risk across a diverse 15 rather than betting on one idea. Once the
measured K_D / k_on / k_off come back, a round-two active-learning pass (e.g. EVOLVEpro) could learn
from those labels to propose the next generation.

# Improving trastuzumab scFv against HER2

Breast Cancer Awareness Hackathon, Zurich, 31 October 2026 (organised by HerCode, challenge by Adaptyv Bio).
Track: **Build a better cancer therapeutic. Trastuzumab scFv optimisation on HER2.**

Source of truth: `FINAL Adaptyv - Zurich hackathon deck (20261030).pdf` (in this folder).

---

## 1. The challenge in one paragraph

Trastuzumab is a main therapy for HER2-overexpressing breast cancer (20-30% of cases). It binds the HER2 extracellular domain close to the membrane (PDB 1N8Z, 2.52 Å). We are given the trastuzumab scFv (242 aa) and must propose **mutated variants that bind HER2 better than the parent (KD 0.82 nM, BLI)**. Each team's best 15 sequences are made and measured by Adaptyv in Lausanne (KD, kon, koff, raw sensorgrams, about 3 weeks later, public on ProteinBase). The track winner gets 45 slots.

## 2. Scope

### In scope
- Sequence-level optimisation of the trastuzumab scFv (point mutations in framework and/or CDRs).
- In-silico ranking of candidates with structure-based confidence metrics (ipSAE from a Boltz-2 cofold with HER2).
- A pool of **15 submitted sequences** plus a **pitch** explaining models used, positions changed and why.

### Out of scope
- De novo binder design (new miniproteins or nanobodies). We are optimising a known binder.
- Wet-lab work (Adaptyv does it), Fc/IgG formatting, ADCC or in-vivo efficacy, clinical claims.
- Changing the scFv architecture (see hard constraints).

### Hard constraints (from the deck, a violation likely means disqualification)
| Rule | Value |
|---|---|
| Mutations from parent | **at least 5** per variant |
| Format | single chain, order **VH - (GGGGS)3 - VL** |
| Length | **up to 260 aa** |
| Submission | **15 variants** per team (45 for the track winner) |

Note: the reverse order (VL-linker-VH) bound at 1.83 nM, so keep VH first as required.

### Judging criteria
1. A creative approach to design optimisation.
2. The tools used and why.
3. **ipSAE of the submitted pool.**
4. Pitch and Q&A.

### Success definition
- **Primary:** a pool that is valid, diverse, and has high ipSAE in Boltz-2 cofolds.
- **Stretch:** at least one variant with KD below 0.82 nM when the lab data arrives.
- Realistic expectation: the parent is already sub-nanomolar, so many variants will be neutral or worse. A portfolio across risk levels is the right strategy (section 5, step 7).

---

## 3. Key data

**Target:** HER2 extracellular domain, 630 aa, Thr23-Thr652, UniProt P04626.
**Structure:** PDB 1N8Z (HER2 ectodomain + trastuzumab Fab). In the scFv, VH = chain B and VL = chain A.

**Parent scFv (242 aa):**

```
VH (120 aa):
EVQLVESGGGLVQPGGSLRLSCAASGFNIKDTYIHWVRQAPGKGLEWVARIYPTNGYTRYADSVKGRFTISADTSKNTAYLQMNSLRAEDTAVYYCSRWGGDGFYAMDYWGQGTLVTVSS
Linker (15 aa):
GGGGSGGGGSGGGGS
VL (107 aa):
DIQMTQSPSSLSASVGDRVTITCRASQDVNTAVAWYQQKPGKAPKLLIYSASFLYSGVPSRFSGSRSGTDFTLTISSLQPEDFATYYCQQHYTTPPTFGQGTKVEIK
```

CDRs (sequence motifs, confirm exact boundaries and numbering in step 2):

| CDR | Motif |
|---|---|
| H1 | GFNIKDTYIH |
| H2 | RIYPTNGYTRYADSVKG |
| H3 | WGGDGFYAMDY |
| L1 | RASQDVNTAVA |
| L2 | SASFLYS |
| L3 | QQHYTTPPT |

---

## 4. Tools: what the deck suggests and how Claude skills map onto it

The deck lists 8 optimisation tools. We use Claude Code with the **Adaptyv protein-design skills** (`adaptyv:*`) to run and orchestrate them (also see `proteinbase.com/protein-design-skills`, 24 skills).

| Deck tool | Role for us | Claude skill available here |
|---|---|---|
| ESM2 / ESM3 | Zero-shot single and double mutant likelihoods (free signal) | `adaptyv:esm` |
| Boltz-2 | Cofold edited scFv + HER2, rank by **ipSAEmin** | `adaptyv:boltz`, `adaptyv:ipsae` |
| Mosaic | Composite objective (affinity + expression + developability) | `adaptyv:mosaic` |
| AbMPNN | Antibody-finetuned ProteinMPNN, backbone-fixed design | `adaptyv:proteinmpnn` (closest; check whether an AbMPNN weight set is usable) |
| AbLang2 | Paired antibody LM, masked scan over framework and CDRs | no dedicated skill, run directly (pip / HF) |
| DSM | Masked-diffusion ESM2, sequence-only edits | no dedicated skill, run directly |
| EVOLVEpro | Few-shot active learning over PLM embeddings | no dedicated skill, only useful if we have measured data (we have almost none, see step 8) |
| ProFam-1 | Family-as-prompt PLM scoring, open weights | no dedicated skill, run directly |

Supporting skills we will use around the models:
- `adaptyv:setup` first-time environment (Modal, biomodals). Run this first.
- `adaptyv:pdb` and `adaptyv:uniprot` fetch 1N8Z and P04626.
- `adaptyv:protein-qc` liabilities (free Cys, N-glycosylation motifs, deamidation, polybasic patches) and biophysical filters.
- `adaptyv:alphafold` / `adaptyv:chai` / `adaptyv:protenix` optional second-opinion cofolds for the final shortlist.
- `adaptyv:campaign-manager` and `adaptyv:protein-design-workflow` to plan, sanity-check pass rates and diagnose failures.
- `adaptyv:binder-design` is for choosing de novo tools, mostly not needed here (we are not designing from scratch).

Other resources from the deck:
- Anthropic's campaign prompt, a good template to adapt: `huggingface.co/datasets/Anthropic/claude-protein-binder-design`, file `prompts/single_target/EGFR.md`.
- Optimising on a small test budget: `adaptyvbio.com/blog/po101` and `po102`.
- ProteinBase (open data, search what already binds a target): `proteinbase.com`.

---

## 5. Step-by-step implementation

Suggested time boxes assume a one-day hackathon. Adjust to the actual schedule.

### Step 0. Environment (about 20 min)
1. Run the `adaptyv:setup` skill: Modal auth, biomodals directory, GPU access.
2. Create a working layout:
   ```
   Protein_Design/
     README.md
     data/        parent.fasta, 1N8Z.pdb, her2_ecd.fasta, scores/*.csv
     scripts/     validate.py, scan_*.py, cofold_*.py, rank.py
     results/     candidates.csv, final_15.csv, cofolds/
     pitch/       slides, notes
   ```
3. Smoke test one Boltz-2 run before anything else so GPU problems surface early.

### Step 1. Validator first (about 20 min)
Write `scripts/validate.py` and use it on every candidate list from now on. It checks:
- Exactly VH + `GGGGSGGGGSGGGGS` + VL composition, VH first, total length at most 260 (parent is 242).
- Hamming distance to the parent at least 5 (count substitutions in aligned positions; if indels are used, define the count conservatively).
- Only the 20 standard amino acids.
- Uniqueness within the pool.

This is the most important script: a rule violation wastes a slot.

### Step 2. Prepare the target and parent (about 30 min)
1. `adaptyv:pdb`: download 1N8Z, extract chain B (VH), chain A (VL) and the HER2 chain.
2. `adaptyv:uniprot`: fetch P04626 and cut Thr23-Thr652 for the ectodomain.
3. Verify the deck sequence matches the PDB chains (note any PDB vs. deck differences).
4. Number the parent sequence (Kabat/IMGT, e.g. with ANARCI) and write out a position table: scFv index, region (FR/CDR), contact with HER2 (distance cutoff, e.g. any heavy atom within 4.5-5 Å), buried or exposed.
5. Decide the cofold target size. Trastuzumab binds the membrane-proximal domain IV, so a truncated HER2 (domain IV plus some context) would make Boltz-2 much faster than the full 630 aa. Choose the cutoff from the 1N8Z contacts and keep it identical across all runs, otherwise ipSAE values are not comparable.

### Step 3. Baseline: score the parent (about 30 min)
1. Cofold the parent scFv + HER2 with Boltz-2 (`adaptyv:boltz`), several seeds.
2. Compute ipSAE (`adaptyv:ipsae`), report ipSAEmin, ipTM, pLDDT.
3. Sanity check: does the predicted pose match 1N8Z (RMSD of the interface, same paratope)? If not, fix the setup (templates, constraints, truncation) before scoring variants. Every later score is relative to this baseline.
4. Also cofold the reverse-order and a shuffled-linker control if time allows, to see how well ipSAE tracks the known 0.82 vs 1.83 nM difference. This is our only calibration signal.

### Step 4. Generate candidate mutations from several independent sources (about 90 min)
Run these in parallel; diversity of evidence is the creative part of the pitch.

| Source | What to do |
|---|---|
| ESM2/ESM3 (`adaptyv:esm`) | Score all 242 x 19 single mutants (about 4.6k), then plausible doubles. Keep mutations whose likelihood is at least as good as the parent residue. |
| AbLang2 | Masked scan over framework and CDRs with the paired model. |
| DSM / ProFam-1 | Independent masked-diffusion / family-prompt proposals. |
| AbMPNN / ProteinMPNN (`adaptyv:proteinmpnn`) | Redesign selected non-core, interface-adjacent positions with the backbone fixed on the Fab-HER2 complex. Keep antigen-contacting key residues fixed. |
| Mosaic (`adaptyv:mosaic`) | Composite loss: Boltz/AF2 interface confidence + stability/expression + developability. Restrict to a small set of positions. |

Guidance for what to mutate:
- Language models mostly propose "more germline-like / more common" residues. That tends to help **stability and expression** and rarely hurts binding, a good low-risk baseline.
- Structure-based design (MPNN, Mosaic) is where **affinity** gains can come from, and also where most failures come from. Limit it to the CDR periphery and second-shell positions.
- Treat CDR-H3 and the central paratope contacts as high risk.
- Do not mutate the Cys residues that form the disulfide bonds, nor the linker (it must remain `(GGGGS)3`).

Output: one tidy table `results/candidates.csv` with columns: variant id, mutations (e.g. `H:S30T`), source tool(s), per-model scores.

### Step 5. Consensus filter (about 30 min)
- Rank mutations by how many independent sources support them (ESM and AbLang2 and structure design agreeing is a strong signal).
- Drop mutations with liabilities using `adaptyv:protein-qc`: new N-glycosylation motifs (N-X-S/T), unpaired Cys, extra deamidation/isomerisation hotspots (NG, DG), large hydrophobic patches, unwanted pI shifts.
- Drop mutations that break the interface in a quick local check (loss of contacts, clashes).

### Step 6. Combine into variants and cofold at scale (about 60-90 min)
1. Build several hundred to a few thousand combinations of at least 5 mutations each, starting from the consensus-ranked list.
2. Cofold with Boltz-2 (seeds fixed, same HER2 construct as the baseline), compute ipSAE with `adaptyv:ipsae`.
3. Keep only variants where ipSAEmin is at or above the parent baseline (and pLDDT is not degraded).
4. Validate every variant with `scripts/validate.py`.
5. Use `adaptyv:campaign-manager` to review the pass rate per stage. If almost nothing passes, loosen thresholds or reduce mutation count per variant, do not add more mutations.

### Step 7. Select the final 15: build a portfolio (about 30 min)
With only 15 slots and a strong parent, diversify by risk instead of picking the top 15 by one score. Suggested split (adjust once we see the data):

| Tier | Slots | Idea |
|---|---|---|
| A. Safe | about 5 | Exactly 5-6 mutations, mostly framework / stability / germline-like, high ESM + AbLang2 agreement. Goal: retain about 1 nM binding. |
| B. Moderate | about 6 | Mix of LM-supported framework changes and 1-3 CDR-periphery affinity candidates from structure design. |
| C. Bold | about 4 | CDR-focused, Mosaic / MPNN-driven, the best shot at KD below 0.82 nM. |

Selection rules: diversity (no two variants share more than about 60-70% of their mutations), every variant passes validation and the QC filters, and the pool-level mean ipSAE is high (the judging metric). Save as `results/final_15.csv`.

### Step 8. Optional: use the small-budget active-learning idea (only if there is data)
EVOLVEpro-style few-shot learning needs measured labels. We only get lab data about 3 weeks after the event, so it is mostly a **future round** (use results on ProteinBase to retrain and design round 2). Check ProteinBase for existing trastuzumab / HER2 measurements first (`proteinbase.com`), they could serve as seed labels.

### Step 9. Pitch (prepare in parallel from step 4 on)
Cover: (1) the problem and why the parent is hard to beat, (2) the tools used and **why each one**, (3) which positions were changed and the evidence per position (table from step 5), (4) the portfolio logic (safe / moderate / bold), (5) ipSAE of the pool vs. the parent baseline, (6) risks and what we would do in round 2. Prepare a one-slide table of the 15 variants with mutations.

---

## 6. Risks and mitigations

| Risk | Mitigation |
|---|---|
| Rule violation (mutation count, length, order) | `validate.py` gate on everything, run it before submitting |
| ipSAE improves but real binding does not (cofold models are imperfect) | Calibrate with the parent and the reverse-order control; use safe-tier variants; do not rely on one predictor (optionally cross-check the top candidates with Chai or Protenix) |
| Improved scores are artefacts of a different HER2 construct | Same construct, seeds and settings for all runs |
| Too many CDR mutations destroy a sub-nM binder | Portfolio approach, keep high-risk variants to a few slots |
| Poor expression / aggregation of the scFv | LM-supported framework changes, `adaptyv:protein-qc` developability checks, avoid new hydrophobic patches |
| Running out of GPU or time | Truncated HER2 construct, small seed counts first, cache results to CSV |

## 7. Deliverables checklist

- [ ] `data/parent.fasta` verified against the deck and 1N8Z
- [ ] `scripts/validate.py` passing on the final pool
- [ ] Baseline Boltz-2 / ipSAE of the parent
- [ ] `results/candidates.csv` with per-tool evidence
- [ ] `results/final_15.csv` (id, full 242-260 aa sequence, mutation list, ipSAEmin, tool provenance)
- [ ] Pitch slides
- [ ] Submission to the organisers in the requested format (confirm the exact format, likely CSV)
- [ ] After the event: lab results from ProteinBase (KD, kon, koff) and a short retrospective

## 8. Open questions to confirm with the organisers
- Exact submission format and deadline.
- How mutations are counted if an insertion/deletion is used (we assume substitutions only).
- Whether the ipSAE used for judging is computed by them (and with which model and HER2 construct) or reported by us.
- Whether the linker or termini (e.g. tags) may be modified within the 260 aa limit (we assume the linker stays `(GGGGS)3`).

## 9. Notes on the Claude-assisted workflow
- Adapt Anthropic's EGFR campaign prompt to HER2/trastuzumab as the master instruction, then let Claude run the skills stage by stage.
- Keep human control on the decisions that matter: which positions are allowed to change, the final 15, and the rule check.
- Log every tool call and its parameters (this is also the material for the "tools used and why" part of the pitch).

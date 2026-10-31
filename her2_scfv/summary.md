# How we estimate K_D: method, research and sources

## 10.1 What is measured and what is predicted
- **Measured:** the parent K_D, 0.82 nM by BLI on Adaptyv's platform. This is the only measured affinity we have.
- **Predicted:** a K_D for each variant, anchored to the parent, reported as a range. These are estimates to be tested in the wet lab, not measurements.

## 10.2 How the method was chosen (literature research)
We searched the literature for computational affinity prediction for protein-protein and antibody-antigen complexes. We kept a method only if it:
1. predicts the change in binding free energy (ddG) for point mutations on a fixed complex,
2. runs on structures we already have (Boltz-2 CIF files),
3. is open or academically licensed, and
4. has been validated against measured binding data on antibodies or protein-protein complexes.

What the literature showed:
- **PRODIGY (absolute affinity from contacts).** A contact-count linear model with Pearson r = 0.73 and RMSE 1.89 kcal/mol across 81 protein-protein complexes (Vangone & Bonvin, eLife 2015; Xue et al., Bioinformatics 2016). Its developers note it was not trained on antibody-antigen complexes, so we use it only as an absolute-scale sanity check.
- **StaB-ddG (deep-learning ddG, ICML 2025).** Transfers folding-stability knowledge from ProteinMPNN to binding. Per-interface Spearman 0.45 on the SKEMPI v2 binding test split. On the authors' split it was not statistically different from FoldX or Flex ddG, and it is over 1,000x faster. Selected as the primary ddG method.
- **FoldX and Rosetta Flex ddG (physics-based).** Flex ddG samples backbone ensembles before scoring and was designed to handle interface mutations, including antibody-antigen interfaces (Barlow et al., 2018). Both need licences, so they are cross-checks when available.
- **Graphinity (equivariant graph neural network, Nature Computational Science 2025).** Test Pearson correlation up to 0.87 on experimental ddG. The authors conclude that current experimental data are insufficient to predict ddG robustly, with orders of magnitude more data likely needed. This is why we widen uncertainty rather than trust any single predictor.
- **Benchmark of eight non-rigorous methods (654 single mutations).** Most performed worse on antibody-antigen complexes than on other complexes. They classify stabilising versus destabilising mutations well but are weaker at quantitative affinity. We therefore use the ddG methods for direction and rough magnitude, not precise numbers.
- **Best antibody-antigen result in that literature:** FoldX combined with molecular dynamics reached r = 0.46, but was the most expensive method.
- **Precedent on this exact system:** an active-learning workflow combined the RDE-Network model with Rosetta Flex ddG on HER2-binding trastuzumab mutants. It found mutants with better binding without any experimental ddG data (arXiv preprint, "Active learning for energy-based antibody optimization and enhanced screening").
- **Boltz-2 affinity head:** trained for small-molecule ligands. It does not predict protein-protein affinity, so we do not use it for K_D.
- **Interface confidence (ipSAE, ipTM):** saturates near 1 for good binders and does not map to K_D. We use it only to discard broken designs.

## 10.3 The estimation method
1. **Anchor.** The parent is the reference. Its free energy is dG = RT ln(K_D). With RT = 0.5925 kcal/mol at 25 °C and K_D = 0.82 nM, dG_parent is about -12.4 kcal/mol.
2. **Predict ddG per variant** with StaB-ddG, using the parent complex as the wild-type structure and each variant complex as the mutant. Add FoldX if licensed.
3. **Sanity check** the absolute scale with PRODIGY on all complexes.
4. **Consensus and uncertainty.** Average the ddG methods. Use an uncertainty of at least 1.0 kcal/mol. Flag any variant whose methods disagree by more than 1.5 kcal/mol.
5. **Convert:** K_D(variant) = 0.82 nM x exp(ddG / RT), with the range from exp(+/- sigma / RT).
6. **Filter.** Exclude variants with low interface confidence or that break the hard rules.

Worked example (illustrative only, not a result): a ddG of +0.5 kcal/mol gives about 1.9 nM, and -0.5 kcal/mol gives about 0.35 nM. One kcal/mol of uncertainty is roughly a 5x range in K_D.

## 10.4 Tools used
| Tool | Source | Role | Status |
|---|---|---|---|
| StaB-ddG | github.com/LDeng0205/StaB-ddG (Deng et al., ICML 2025) | primary ddG | to install and verify |
| PRODIGY | pip prodigy-prot (Vangone & Bonvin 2015; Xue et al. 2016) | absolute dG sanity check | to install and verify |
| FoldX | foldxsuite (licence required) | ddG cross-check | optional |
| Rosetta Flex ddG | Barlow et al. 2018 (licence required) | ddG cross-check | skipped for now |
| ipSAE | Dunbrack, bioRxiv 2025.02.10.637595 | interface-confidence filter | done (results.md) |
| ESM2, AbLang2 | earlier sequence scans | plausibility prior, flags only | done |
| Boltz-2 | modal_boltz.py | structure of each complex | done (cofold/) |

## 10.5 Limits we state openly
- Antibody-antigen ddG predictors reach only modest accuracy (r about 0.4 to 0.5), with roughly 1 kcal/mol error. That is about a 5x range in K_D.
- Variants carry 5 to 7 simultaneous mutations. Effects are not fully additive, so their ranges are wider.
- Structures are Boltz predictions, not crystal structures. Errors in the complex propagate into the ddG.
- The ddG methods were trained on limited experimental data, so they can overfit.
- The measured K_D from Adaptyv is the only decisive answer.

## 10.6 References
- Vangone A., Bonvin A.M.J.J. (2015). Contacts-based prediction of binding affinity in protein-protein complexes. eLife. doi:10.7554/eLife.07454
- Xue L.C. et al. (2016). PRODIGY: a web server for predicting the binding affinity of protein-protein complexes. Bioinformatics. doi:10.1093/bioinformatics/btw514
- Deng A. et al. (2025). Predicting mutational effects on protein binding from folding energy. ICML 2025, PMLR 267:13129-13151. arXiv:2507.05502. Code: github.com/LDeng0205/StaB-ddG
- Barlow K.A. et al. (2018). Flex ddG: Rosetta ensemble-based estimation of changes in protein-protein binding affinity upon mutation. PMC5980710 (bioRxiv 221689)
- Graphinity (2025). Investigating the volume and diversity of data needed for generalizable antibody-antigen ddG prediction. Nature Computational Science 5, 635-647. doi:10.1038/s43588-025-00823-8
- Active learning for energy-based antibody optimization and enhanced screening. arXiv preprint (HER2 trastuzumab case study). Author list to be confirmed before citing.
- Dunbrack R.L. (2025). Res ipSAE loquuntur: what's wrong with AlphaFold's ipTM score and how to fix it. bioRxiv doi:10.1101/2025.02.10.637595
- Benchmark of eight non-rigorous methods on 654 single mutations. Journal and author list to be confirmed before citing.

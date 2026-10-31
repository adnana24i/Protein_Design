#!/usr/bin/env bash
# Builds the two Boltz co-fold input files (parent + AH06/D31E variant), each = scFv chain + HER2 ECD chain.
# Run this ONCE from your real VS Code terminal (it fetches HER2 from UniProt; needs internet).
set -euo pipefail
cd "$(dirname "$0")"

# scFv sequences (VH - (GGGGS)3 - VL), 242 aa each
PARENT="EVQLVESGGGLVQPGGSLRLSCAASGFNIKDTYIHWVRQAPGKGLEWVARIYPTNGYTRYADSVKGRFTISADTSKNTAYLQMNSLRAEDTAVYYCSRWGGDGFYAMDYWGQGTLVTVSSGGGGSGGGGSGGGGSDIQMTQSPSSLSASVGDRVTITCRASQDVNTAVAWYQQKPGKAPKLLIYSASFLYSGVPSRFSGSRSGTDFTLTISSLQPEDFATYYCQQHYTTPPTFGQGTKVEIK"
# Variant = AH06 (D102W,M107F,D108A,Y109L in CDR-H3) + D31E (CDR-H1) = 5 mutations from parent
VARIANT="EVQLVESGGGLVQPGGSLRLSCAASGFNIKETYIHWVRQAPGKGLEWVARIYPTNGYTRYADSVKGRFTISADTSKNTAYLQMNSLRAEDTAVYYCSRWGGWGFYAFALWGQGTLVTVSSGGGGSGGGGSGGGGSDIQMTQSPSSLSASVGDRVTITCRASQDVNTAVAWYQQKPGKAPKLLIYSASFLYSGVPSRFSGSRSGTDFTLTISSLQPEDFATYYCQQHYTTPPTFGQGTKVEIK"

echo "Fetching HER2 (UniProt P04626)..."
curl -sS "https://rest.uniprot.org/uniprotkb/P04626.fasta" -o her2_full.fasta
HER2=$(python3 -c "s=''.join(l.strip() for l in open('her2_full.fasta') if not l.startswith('>')); print(s[22:652])")

echo "scFv parent length:  ${#PARENT}  (expect 242)"
echo "scFv variant length: ${#VARIANT}  (expect 242)"
echo "HER2 ECD length:     ${#HER2}  (expect 630)"

printf ">scfv\n%s\n>her2\n%s\n" "$PARENT" "$HER2" > parent_complex.fasta
printf ">scfv\n%s\n>her2\n%s\n" "$VARIANT" "$HER2" > ah06_d31e_complex.fasta
echo "Wrote parent_complex.fasta and ah06_d31e_complex.fasta"

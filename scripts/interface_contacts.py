"""Per-residue distance to HER2 for the trastuzumab VH/VL in PDB 1N8Z.

Usage: python3 -I scripts/interface_contacts.py data/pdb1n8z.ent data/interface_contacts.csv
Numbering: VH/VL use PDB (sequential) numbering; scFv index = VH n, or 135 + VL n.
CDRs are located by sequence motif (see CDRS), not by a formal numbering scheme.
"""
import csv, math, sys

AA3 = {'ALA':'A','ARG':'R','ASN':'N','ASP':'D','CYS':'C','GLN':'Q','GLU':'E','GLY':'G','HIS':'H','ILE':'I',
       'LEU':'L','LYS':'K','MET':'M','PHE':'F','PRO':'P','SER':'S','THR':'T','TRP':'W','TYR':'Y','VAL':'V'}
CUTOFF = 4.5
DOMAINS = {'B': ('VH', 120, 0), 'A': ('VL', 107, 135)}  # chain -> (name, length, scFv offset)
CDRS = {'H1': 'GFNIKDTYIH', 'H2': 'RIYPTNGYTRYADSVKG', 'H3': 'WGGDGFYAMDY',
        'L1': 'RASQDVNTAVA', 'L2': 'SASFLYS', 'L3': 'QQHYTTPPT'}

pdb, out = sys.argv[1], sys.argv[2]
res = {}
for line in open(pdb):
    if line.startswith('ATOM'):
        key = (line[21], int(line[22:26]), AA3[line[17:20].strip()])
        res.setdefault(key, []).append((float(line[30:38]), float(line[38:46]), float(line[46:54])))

her2 = [(k, xyz) for k, xyz in res.items() if k[0] == 'C']
seq = {ch: ''.join(k[2] for k in sorted(res) if k[0] == ch and k[1] <= n) for ch, (_, n, _) in DOMAINS.items()}
cdr_of = {}
for name, motif in CDRS.items():
    ch = 'B' if name[0] == 'H' else 'A'
    start = seq[ch].index(motif) + 1
    for i in range(start, start + len(motif)):
        cdr_of[(ch, i)] = name

rows = []
for (ch, n, aa), xyz in sorted(res.items()):
    if ch not in DOMAINS or n > DOMAINS[ch][1]:
        continue
    best, partner = 99.0, None
    for (_, hn, haa), hxyz in her2:
        d = min(math.dist(a, b) for a in xyz for b in hxyz)
        if d < best:
            best, partner = d, hn
    name, _, off = DOMAINS[ch]
    rows.append({'domain': name, 'pdb_chain': ch, 'pdb_resnum': n, 'aa': aa, 'scfv_index': off + n,
                 'region': cdr_of.get((ch, n), 'FR'), 'min_dist_to_her2': round(best, 2),
                 'nearest_her2_pdb_resnum': partner, 'nearest_her2_uniprot': partner + 22,
                 'contact_4.5A': best <= CUTOFF})

with open(out, 'w', newline='') as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0]))
    w.writeheader()
    w.writerows(rows)
print(len(rows), 'residues;', sum(r['contact_4.5A'] for r in rows), 'contacts ->', out)

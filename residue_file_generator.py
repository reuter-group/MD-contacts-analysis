#!/usr/bin/env python

import sys 
import MDAnalysis


SEGID= 'PROA'

if len(sys.argv) < 2:
    print("Usage: ./residue_file_generator.py PSF_FILENAME [options]")
    print("  options:")
    print("  -a Only aromatic residues, i.e. PHE, TYR, TRP")
    sys.exit()

u = MDAnalysis.Universe(sys.argv[1])
atoms = u.select_atoms("segid {}".format(SEGID))

with open(sys.argv[1]+'_candidates.txt', 'w') as f:
    if len(sys.argv) == 3 and sys.argv[2] == '-a':
        for r in atoms.residues:
            if r.resname.upper() in ['PHE', 'TYR', 'TRP']:
                f.write("{}{}\n".format(r.resname.upper(), r.resid))
    else:
        for r in atoms.residues:
            f.write("{}{}\n".format(r.resname, r.resid))
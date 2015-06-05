#!/usr/bin/env python
#
# Interface with iPython notebook.
#
# Tristan BEREAU (2015)

from __future__ import print_function

import chemview
from chemview import MolecularViewer, enable_notebook, RepresentationViewer
enable_notebook()
import numpy as np

def draw_mol(mol):
    num_atoms = mol.GetConformer().GetNumAtoms()
    coordinates = []
    atom_types = []
    bonds = []
    for i in xrange(num_atoms):
        coordinates.append([j for j in mol.GetConformer().GetAtomPosition(i)])
        atom_types.append(mol.GetAtomWithIdx(i).GetSymbol())
    for j in xrange(len(mol.GetBonds())):
        bond_j = mol.GetBondWithIdx(j)
        bonds.append((bond_j.GetBeginAtomIdx(),bond_j.GetEndAtomIdx()))
    mv2 = MolecularViewer(np.array(coordinates),
                         topology={'atom_types': atom_types,
                                   'bonds': bonds})
    mv2.points()
    #mv2.lines()
    mv2.ball_and_sticks()
    return mv2

def draw_cg(gro, itp):
    coordinates = []
    sizes = []
    colors = []
    atom_types = []
    color_mapping = {
        'Qda': 0x3333CC,
        'Qa' : 0x0000CC,
        'Qd' : 0x0066FF,
        'Q0' : 0x00CCFF,
        'P5' : 0xFF3333,
        'P4' : 0xCC3399,
        'P3' : 0xFF3399,
        'P2' : 0xFF6600,
        'P1' : 0xFF3300,
        'Nda': 0xCC6600,
        'Na' : 0xFF9933,
        'Nd' : 0xFFFF00,
        'N0' : 0xFFFF66,
        'C5' : 0x333333,
        'C4' : 0x666666,
        'C3' : 0x999999,
        'C2' : 0xCCCCCC,
        'C1' : 0xFFFFFF,
    }
    try:
        f = open(gro,'r')
        s = f.readlines()
        f.close()
    except IOError as e:
        print("Error. Can't read file",gro)
        print(e)
        exit(1)
    num_atoms = int(s[1])
    for i in xrange(2,2+num_atoms):
        line = s[i].split()
        coordinates.append([
                float(line[3]),float(line[4]),float(line[5])
            ])
    try:
        f = open(itp,'r')
        s = f.readlines()
        f.close()
    except IOError as e:
        print("Error. Can't read file",itp)
        print(e)
        exit(1)
    parse_on = False
    parse_off = False
    for i in xrange(len(s)):
        if parse_on and not parse_off:
            if s[i][0] != ";":
                line = s[i].split()
                if len(line) > 1:
                    atom_type = line[1]
                    satom_type = atom_type
                    if atom_type[0] == "S":
                        sizes.append(0.14)
                        atom_type = atom_type[1:]
                    else:
                        sizes.append(0.18)
                    if atom_type in color_mapping:
                        colors.append(color_mapping[atom_type])
                    else:
                        print("Error. Unknown atom type",atom_type)
                    print(satom_type)
        if "[atoms]" in s[i]:
            parse_on = True
        if "[constraints]" in s[i] or "[bonds]" in s[i]:
            parse_off = True
    coordinates = np.array(coordinates, 'float32')
    rv = RepresentationViewer()
    spheres_id  = rv.add_representation('spheres', {'coordinates': coordinates,
                                                    'colors': colors,
                                                    'radii': sizes,
                                                    'resolution': 32})
    return rv

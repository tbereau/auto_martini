#!/usr/bin/env python
#
# Automatic MARTINI mapping and parametrization of small organic molecules.
#
# Tristan BEREAU (2014)

from __future__ import print_function

import sys
import os
import argparse
import math
import numpy as np
import six
from numpy import arctan2
from collections import Counter
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import ChemicalFeatures
from rdkit.Chem import rdchem
from rdkit.Chem import rdmolops
from rdkit.Chem import rdmolfiles
from rdkit.Chem import rdMolTransforms
from rdkit.Chem import rdMolDescriptors
from rdkit import RDConfig
import requests
import BeautifulSoup
from itertools import chain
from collections import defaultdict
from operator import itemgetter


# For feature extraction
fdefName = os.path.join(RDConfig.RDDataDir, 'BaseFeatures.fdef')
factory = ChemicalFeatures.BuildFeatureFactory(fdefName)


# Measured octanol/water free energies from MARTINI.
# Data used TI, not partitioning of Marrink et al. JPCB 2007.
delta_f_types = {
    'Qda': -15.04,
    'Qa': -15.04,
    'Qd': -15.04,
    'Q0': -22.35,
    'P5': -8.88,
    'P4': -9.30,
    'P3': -8.81,
    'P2': -3.85,
    'P1': -2.26,
    'Nda': 2.49,
    'Na': 2.49,
    'Nd': 2.49,
    'N0': 4.22,
    'C5': 6.93,
    'C4': 10.14,
    'C3': 12.26,
    'C2': 13.74,
    'C1': 14.20,
}


# Parameters
# CG Bead vdw radius (in Angstroem)
rvdw = 4.7 / 2.
rvdwAromatic = 4.3 / 2.
rvdwCross = 0.5 * (rvdw + rvdwAromatic)
# Optimized parameters
offsetBeadWeight = 50
offsetBeadAromaticWeight = 20.
lonelyAtomPenalize = 0.20
bdBdOverlapCoeff = 9.0
atInBdCoeff = 0.9


def gen_molecule_smi(smi):
    """Generate mol object from smiles string"""
    molecule = Chem.MolFromSmiles(smi)
    molecule = Chem.AddHs(molecule)
    AllChem.EmbedMolecule(molecule, useRandomCoords=True)
    try:
        AllChem.UFFOptimizeMolecule(molecule)
    except ValueError as e:
        print(e)
        exit(1)
    return molecule


def gen_molecule_sdf(sdf):
    """Generate mol object from SD file"""
    suppl = Chem.SDMolSupplier(sdf)
    if len(suppl) > 1:
        print("Error. Only one molecule may be provided.")
        exit(1)
    molecule = ''
    for molecule in suppl:
        if molecule is None:
            print("Error. Can't read molecule.")
            exit(1)
    return molecule


def print_header(arguments):
    """Print topology header"""
    print(";;;; GENERATED WITH auto-martini")
    if arguments.smi:
        print("; INPUT SMILES:", arguments.smi)
    else:
        print(";", arguments.sdf)
    print("; Tristan Bereau (2014)")
    print("")
    print('[moleculetype]')
    print('; molname       nrexcl')
    print('  {:5s}         2'.format(arguments.molname))
    print('')
    print('[atoms]')
    print('; id    type    resnr   residu  atom    cgnr    charge  smiles')
    return


def letter_occurrences(string):
    """Count letter occurence"""
    frequencies = defaultdict(lambda: 0)
    for character in string:
        if character.isalnum():
            frequencies[character.upper()] += 1
    return frequencies


def get_charge(molecule):
    """Get net charge of molecule"""
    return Chem.rdmolops.GetFormalCharge(molecule)


def get_hbond_a(features):
    """Get Hbond acceptor information"""
    hbond = []
    for feat in features:
        if feat.GetFamily() == "Acceptor":
            for i in feat.GetAtomIds():
                if i not in hbond:
                    hbond.append(i)
    return hbond


def get_hbond_d(features):
    """Get Hbond donor information"""
    hbond = []
    for feat in features:
        if feat.GetFamily() == "Donor":
            for i in feat.GetAtomIds():
                if i not in hbond:
                    hbond.append(i)
    return hbond


def extract_features(molecule):
    """Extract features of mol"""
    features = factory.GetFeaturesForMol(molecule)
    # if len(feats) == 0:
    # print "Error. Can't extract molecular features."
    # exit(1)
    return features


def output_xyz(arguments, molecule):
    """Ouput XYZ file of molecule object"""
    num_atoms = molecule.GetConformer().GetNumAtoms()
    if arguments.xyz[-4:] != ".xyz":
        arguments.xyz += ".xyz"
    try:
        with open(arguments.xyz, 'w') as f:
            f.write(str(num_atoms) + "\n")
            f.write(" " + arguments.xyz[:-4] + "\n")
            for i in range(num_atoms):
                f.write("{:2s}    {:7.4f} {:7.4f} {:7.4f}\n".format(
                    molecule.GetAtomWithIdx(i).GetSymbol(),
                    molecule.GetConformer().GetAtomPosition(i)[0],
                    molecule.GetConformer().GetAtomPosition(i)[1],
                    molecule.GetConformer().GetAtomPosition(i)[2]))
            f.write("\n")
            f.close()
    except IOError:
        print("Can't write to file " + arguments.xyz)
        exit(1)
    return


def output_gro(arguments, beads, bead_names):
    """Output GRO file of CG structure"""
    num_beads = len(beads)
    if len(beads) != len(bead_names):
        print("Error. Incompatible number of beads and bead names.")
        exit(1)
    if arguments.gro[-4:] != ".gro":
        arguments.gro += ".gro"
    try:
        with open(arguments.gro, 'w') as f:
            f.write("{:s} generated from {:s}\n".format(
                arguments.molname, os.path.basename(__file__)))
            f.write("{:5d}\n".format(num_beads))
            for i in range(num_beads):
                f.write("{:5d}{:<6s} {:3s}{:5d}{:8.3f}{:8.3f}{:8.3f}\n".format(
                    i + 1, arguments.molname, bead_names[i], i + 1, beads[i][0] / 10.,
                    beads[i][1] / 10., beads[i][2] / 10.))
            f.write("{:10.5f}{:10.5f}{:10.5f}\n".format(10., 10., 10.))
            f.close()
    except IOError:
        print("Can't write to file " + arguments.gro)
        exit(1)
    return


def output_pdb(molecule, cgbeads):
    """Output PDB file of AA/CG structure"""
    for i in range(len(cgbeads)):
        ati = molecule.GetAtomWithIdx(i)
        pdbres = Chem.rdchem.AtomPDBResidueInfo(ati)
        print(pdbres.GetResidueName())
        print(ati.GetSmarts())
    pw = Chem.rdmolfiles.PDBWriter(args.pdb)
    Chem.rdmolfiles.PDBWriter.write(pw, molecule)
    pw.close()
    return


def get_ring_atoms(features):
    """Get ring atoms"""
    ringatoms = []
    for feat in features:
        if feat.GetType() in ["RH6_6", "RH5_5", "RH4_4", "RH3_3",
                              "Arom5", "Arom6", "Arom7", "Arom8"]:
            new_ring = []
            for at in feat.GetAtomIds():
                new_ring.append(at)
            if new_ring not in ringatoms:
                ringatoms.append(new_ring)
    if args.verbose:
        print("; ring atoms:", ringatoms)
    return ringatoms


def gaussian_overlap(molecule, bead1, bead2, ringatoms):
    """"Returns overlap coefficient between two gaussians
    given distance dist"""
    conf = molecule.GetConformer()
    dist = Chem.rdMolTransforms.GetBondLength(conf, bead1, bead2)
    sigma = rvdw
    if bead1 in ringatoms and bead2 in ringatoms:
        sigma = rvdwAromatic
    if bead1 in ringatoms and bead2 not in ringatoms or \
       bead1 not in ringatoms and bead2 in ringatoms:
        sigma = rvdwCross
    return bdBdOverlapCoeff * math.exp(-dist ** 2 / 4. / sigma ** 2)


def atoms_in_gaussian(molecule, bead_id, ringatoms):
    """Returns weighted sum of atoms contained in bead bead_id"""
    weight_sum = 0.0
    conf = molecule.GetConformer()
    sigma = rvdw
    lumped_atoms = []
    if bead_id in ringatoms:
        sigma = rvdwAromatic
    for i in range(conf.GetNumAtoms()):
        dist_bd_at = Chem.rdMolTransforms.GetBondLength(conf, i, bead_id)
        if dist_bd_at < sigma:
            lumped_atoms.append(i)
        weight_sum -= molecule.GetAtomWithIdx(i).GetMass() * math.exp(-dist_bd_at ** 2 / 2 / sigma ** 2)
    return atInBdCoeff * weight_sum, lumped_atoms


def penalize_lonely_atoms(molecule, lumped_atoms):
    """Penalizes configuration if atoms aren't included
    in any CG bead"""
    weight_sum = 0.0
    conf = molecule.GetConformer()
    for i in range(conf.GetNumAtoms()):
        if i not in lumped_atoms:
            weight_sum += molecule.GetAtomWithIdx(i).GetMass()
    return lonelyAtomPenalize * weight_sum


def eval_gaussian_interac(molecule, list_beads, ringatoms):
    """From collection of CG beads placed on mol, evaluate
    objective function of interacting beads"""
    weight_sum = 0.0
    # Offset energy for every new CG bead.
    # Distinguish between aromatics and others.
    num_aromatics = 0
    lumped_atoms = []
    for i in list_beads:
        if i in ringatoms:
            num_aromatics += 1
    weight_sum += offsetBeadWeight * (len(list_beads) - num_aromatics) + \
        offsetBeadAromaticWeight * num_aromatics
    # Repulsive overlap between CG beads
    for i in range(len(list_beads)):
        for j in range(i + 1, len(list_beads)):
            weight_sum += gaussian_overlap(molecule, list_beads[i], list_beads[j], ringatoms)
    # Attraction between atoms nearby to CG bead
    for i in range(len(list_beads)):
        weight, lumped = atoms_in_gaussian(molecule, list_beads[i], ringatoms)
        weight_sum += weight
        for j in lumped:
            if j not in lumped_atoms:
                lumped_atoms.append(j)
    # Penalty for excluding atoms
    weight_sum += penalize_lonely_atoms(molecule, lumped_atoms)
    return weight_sum


def enumerate_seq(heavy_atoms, depth):
    """Enumerate all sequences of length depth among heavy atoms"""
    seq = [[heavy_atoms[0]] * depth]
    while seq[-1] != [heavy_atoms[-1]] * depth:
        last_ele = seq[-1]
        new_ele = []
        last_index = 1
        toggle_next = True
        while last_index <= depth:
            if toggle_next:
                if last_ele[-last_index] == heavy_atoms[-1]:
                    # find left-most number that's not == heavy_atoms[-1]
                    lix = 0
                    for lix in range(1, depth + 1):
                        if last_ele[-lix] != heavy_atoms[-1]:
                            break
                    for lx in range(0, lix - 1):
                        new_ele = [heavy_atoms[heavy_atoms.index(last_ele[-lix])+1]] + new_ele
                    last_index = lix - 1
                    toggle_next = True
                else:
                    new_ele = [heavy_atoms[heavy_atoms.index(
                        last_ele[-last_index]) + 1]] + new_ele
                    toggle_next = False
            else:
                new_ele = [last_ele[-last_index]] + new_ele
            last_index += 1
        if not toggle_next:
            seq.append(new_ele)
    return seq


def find_bead_pos(molecule, ringatoms):
    """Try out all possible combinations of CG beads
    up to threshold number of beads per atom. find
    arrangement with best energy score. Return all
    possible arrangements sorted by energy score."""
    # number of atoms in mol
    conf = molecule.GetConformer()
    num_atoms = conf.GetNumAtoms()
    # List of heavy atoms
    list_heavy_atoms = []
    print("-----------------------")
    if args.verbose:
        print("; Heavy atoms:")
        print("; ", end="")
    for i in range(num_atoms):
        if molecule.GetAtomWithIdx(i).GetSymbol() != "H":
            list_heavy_atoms.append(i)
            if args.verbose:
                print(molecule.GetAtomWithIdx(i).GetSymbol(), end="")
    if args.verbose:
        print("")
    if len(list_heavy_atoms) == 0:
        print("Error. No heavy atom found.")
        exit(1)
    if len(list_heavy_atoms) == 1:
        # Put one CG bead on the one heavy atom.
        best_trial_comb = enumerate_seq(list_heavy_atoms, 1)[0]
        avg_pos = [[conf.GetAtomPosition(best_trial_comb[0])[j] for j in range(3)]]
        return best_trial_comb, avg_pos
    if len(list_heavy_atoms) > 25:
        print("Error. Exhaustive enumeration can't handle large molecules.")
        print(("Number of heavy atoms:", len(list_heavy_atoms)))
        exit(1)
    ringatoms_flat = list(chain.from_iterable(ringatoms))
    # List of bonds between heavy atoms
    list_bonds = []
    for i in range(len(list_heavy_atoms)):
        for j in range(i + 1, len(list_heavy_atoms)):
            if molecule.GetBondBetweenAtoms(list_heavy_atoms[i], list_heavy_atoms[j]) is not None:
                list_bonds.append([list_heavy_atoms[i], list_heavy_atoms[j]])
    # Max number of beads. At most 2.5 heavy atoms per bead.
    max_beads = int(len(list_heavy_atoms) / 2.)
    # Collect all possible combinations of bead positions
    best_trial_comb = []
    list_trial_comb = []
    ene_best_trial = 1e6
    last_best_trial_comb = []

    # Keep track of all combinations and scores
    list_combs = []
    list_energies = []

    # Heavy atom coordinates
    heavyatom_coords = get_heavy_atom_coords(molecule)

    for num_beads in range(1, max_beads + 1):
        # Use recursive function to loop through all possible
        # combinations of CG bead positions.
        seq_one_beads = enumerate_seq(list_heavy_atoms, num_beads)

        combs = []
        energies = []

        # Trial positions: any heavy atom
        for i in range(len(seq_one_beads)):
            trial_comb = seq_one_beads[i]
            # Check for beads at the same place
            count = Counter(trial_comb)
            all_different = True
            for val in count.values():
                if val != 1:
                    all_different = False
                    break
            if all_different:
                # Check for beads linked by chemical bond (except in rings)
                acceptable_trial = True
                bonds_in_rings = [0] * len(ringatoms)
                for bi in range(len(trial_comb)):
                    for bj in range(bi + 1, len(trial_comb)):
                        if [trial_comb[bi], trial_comb[bj]] in list_bonds \
                                or [trial_comb[bj], trial_comb[bi]] in list_bonds:
                            bond_in_ring = False
                            for r in range(len(ringatoms)):
                                if trial_comb[bi] in ringatoms[r] and trial_comb[bj] in ringatoms[r]:
                                    bonds_in_rings[r] += 1
                                    bond_in_ring = True
                            if not bond_in_ring:
                                acceptable_trial = False
                                break
                if acceptable_trial:
                    # Don't allow bonds between atoms of the same ring.
                    for bir in range(len(bonds_in_rings)):
                        if bonds_in_rings[bir] > 0:
                            acceptable_trial = False
                if acceptable_trial:
                    # Check for two terminal beads linked by only one atom
                    for bi in range(len(trial_comb)):
                        for bj in range(bi + 1, len(trial_comb)):
                            if ([item for sublist in list_bonds for item in
                                 sublist].count(trial_comb[bi]) == 1) and ([item for sublist
                                                                           in list_bonds for item in sublist].count(
                                    trial_comb[bj]) == 1):
                                # Both beads are on terminal atoms. Block contribution
                                # if the two terminal atoms are linked to the same atom. 
                                partneri = ''
                                partnerj = ''
                                for bond in list_bonds:
                                    if bond[0] == trial_comb[bi]:
                                        partneri = bond[1]
                                    if bond[1] == trial_comb[bi]:
                                        partneri = bond[0]
                                    if bond[0] == trial_comb[bj]:
                                        partnerj = bond[1]
                                    if bond[1] == trial_comb[bj]:
                                        partnerj = bond[0]
                                if partneri == partnerj:
                                    acceptable_trial = False
                                    break
                if acceptable_trial:
                    # Do the energy evaluation
                    trial_ene = eval_gaussian_interac(molecule, trial_comb, ringatoms_flat)
                    combs.append(trial_comb)
                    energies.append(trial_ene)
                    if args.verbose:
                        print(";", trial_comb, trial_ene)
                    # Make sure all atoms within one bead would be connected
                    if all_atoms_in_beads_connected(trial_comb,
                       heavyatom_coords, list_heavy_atoms, list_bonds):
                        # Accept the move
                        if trial_ene < ene_best_trial:
                            ene_best_trial = trial_ene
                            best_trial_comb = sorted(trial_comb)
                        # Get bead positions
                        beadpos = [[0]*3 for l in range(len(trial_comb))]
                        for l in range(len(trial_comb)):
                            beadpos[l] = [conf.GetAtomPosition(sorted(trial_comb)[l])[m] for m in range(3)]
                        # Store configuration
                        list_trial_comb.append([trial_comb, beadpos, trial_ene])
        # print best_trial_comb
        if last_best_trial_comb == best_trial_comb:
            break
        last_best_trial_comb = best_trial_comb
        list_combs.append(combs)
        list_energies.append(energies)
    if args.verbose:
        for at in best_trial_comb:
            print("; CG bead:", at)
        print("; with energy:", ene_best_trial)
    sorted_combs = np.array(sorted(list_trial_comb, key=itemgetter(2)))
    return sorted_combs[:, 0], sorted_combs[:, 1]


def get_heavy_atom_coords(molecule):
    """Extract atomic coordinates of heavy atoms in molecule mol"""
    heavyatom_coords = []
    conf = molecule.GetConformer()
    # number of atoms in mol
    num_atoms = molecule.GetConformer().GetNumAtoms()
    for i in range(num_atoms):
        if molecule.GetAtomWithIdx(i).GetSymbol() != "H":
            heavyatom_coords.append(np.array(
                [conf.GetAtomPosition(i)[j] for j in range(3)]))
    return heavyatom_coords


def get_cg_bead_coords(molecule, cgbeads, avg_pos, ringatoms_flat):
    """Extract coordinates of CG beads"""
    # CG beads are averaged over best trial combinations for all
    # non-aromatic atoms.
    cgbead_coords = []
    conf = molecule.GetConformer()
    for i in range(len(cgbeads)):
        if cgbeads[i] in ringatoms_flat:
            cgbead_coords.append(np.array([conf.GetAtomPosition(cgbeads[i])[j] for j in range(3)]))
        else:
            # Use average
            cgbead_coords.append(np.array(avg_pos[i]))
    return cgbead_coords


def all_atoms_in_beads_connected(trial_comb, heavyatom_coords, list_heavy_atoms, bondlist):
    """Make sure all atoms within one CG bead are connected to at least
    one other atom in that bead"""
    # Bead coordinates are given by heavy atoms themselves
    cgbead_coords = []
    for i in range(len(trial_comb)):
        cgbead_coords.append(heavyatom_coords[list_heavy_atoms.index(trial_comb[i])])
    voronoi = voronoi_atoms(cgbead_coords, heavyatom_coords)
    for i in range(len(trial_comb)):
        cg_bead = trial_comb[i]
        num_atoms = voronoi.values().count(voronoi[list_heavy_atoms.index(cg_bead)])
        # sub-part of bond list that only contains atoms within CG bead
        sub_bond_list = []
        for j in range(len(bondlist)):
            if voronoi[list_heavy_atoms.index(bondlist[j][0])] == voronoi[list_heavy_atoms.index(cg_bead)] and \
               voronoi[list_heavy_atoms.index(bondlist[j][1])] == voronoi[list_heavy_atoms.index(cg_bead)]:
                sub_bond_list.append(bondlist[j])
        num_bonds = len(sub_bond_list)
        if num_bonds < num_atoms - 1 or num_atoms == 1:
            return False
    return True


def voronoi_atoms(cgbead_coords, heavyatom_coords):
    """Partition all atoms between CG beads"""
    partitioning = {}
    for j in range(len(heavyatom_coords)):
        if j not in partitioning.keys():
            # Voronoi to check whether atom is closest to bead
            bead_at = -1
            dist_bead_at = 1000
            for k in range(len(cgbead_coords)):
                distk = np.linalg.norm(cgbead_coords[k] - heavyatom_coords[j])
                if distk < dist_bead_at:
                    dist_bead_at = distk
                    bead_at = k
            partitioning[j] = bead_at
    if len(cgbead_coords) > 1:
        # Book-keeping of closest atoms to every bead
        closest_atoms = {}
        for i in range(len(cgbead_coords)):
            closest_atom = -1
            closest_dist = 10000.0
            for j in range(len(heavyatom_coords)):
                dist_bead_at = np.linalg.norm(cgbead_coords[i] - heavyatom_coords[j])
                if dist_bead_at < closest_dist:
                    closest_dist = dist_bead_at
                    closest_atom = j
            if closest_atom == -1:
                print("Error. Can't find closest atom to bead", i)
                exit(1)
            closest_atoms[i] = closest_atom
        # If one bead has only one heavy atom, include one more
        for i in partitioning.values():
            if sum(x == i for x in partitioning.values()) == 1:
                # Find bead
                lonely_bead = i
                # Voronoi to find closest atom
                closest_bead = -1
                closest_bead_dist = 10000.0
                for j in range(len(heavyatom_coords)):
                    if partitioning[j] != lonely_bead:
                        dist_bead_at = np.linalg.norm(cgbead_coords[lonely_bead] - heavyatom_coords[j])
                        # Only consider if it's closer, not a CG bead itself, and
                        # the CG bead it belongs to has more than one other atom.
                        if dist_bead_at < closest_bead_dist and \
                           j != closest_atoms[partitioning[j]] and \
                           sum(x == partitioning[j] for x in partitioning.values()) > 2:
                            closest_bead = j
                            closest_bead_dist = dist_bead_at
                if closest_bead == -1:
                    print("Error. Can't find an atom close to atom", lonely_bead)
                    exit(1)
                partitioning[closest_bead] = lonely_bead
    return partitioning


def substruct2smi(molecule, partitioning, cg_bead, cgbeads, ringatoms):
    """Substructure to smiles conversion; also output Wildman-Crippen log_p;
         and charge of group."""
    frag = rdchem.EditableMol(molecule)
    num_atoms = molecule.GetConformer().GetNumAtoms()
    # First delete all hydrogens
    for i in range(num_atoms):
        if molecule.GetAtomWithIdx(i).GetSymbol() == "H":
            # find atom from coordinates
            submol = frag.GetMol()
            for j in range(submol.GetConformer().GetNumAtoms()):
                if molecule.GetConformer().GetAtomPosition(i)[0] == \
                        submol.GetConformer().GetAtomPosition(j)[0]:
                    frag.RemoveAtom(j)
    # Identify atoms involved in same ring as cg_bead (only one ring)
    atoms_in_ring = []
    for ring in ringatoms:
        if cgbeads[cg_bead] in ring:
            atoms_in_ring = ring[:]  # CHANGED
            break
    # Add atoms off the ring that belong to the fragment.
    for atom in atoms_in_ring:
        # if atom in cg_beads:
        if atom == cgbeads:
            for atp in partitioning.keys():
                if partitioning[atp] == partitioning[atom] and atp not in atoms_in_ring:
                    atoms_in_ring.append(atp)
    # Then heavy atoms that aren't part of the CG bead (except those
    # involved in the same ring).
    for i in partitioning.keys():
        if partitioning[i] != cg_bead and i not in atoms_in_ring:
            # find atom from coordinates
            submol = frag.GetMol()
            for j in range(submol.GetConformer().GetNumAtoms()):
                if molecule.GetConformer().GetAtomPosition(i)[0] == \
                        submol.GetConformer().GetAtomPosition(j)[0]:
                    frag.RemoveAtom(j)
    # Wildman-Crippen log_p
    wc_log_p = rdMolDescriptors.CalcCrippenDescriptors(frag.GetMol())[0]
    tmpfile = 'tmp-auto-martini.smi'
    sw = Chem.rdmolfiles.SmilesWriter(tmpfile, nameHeader='')
    # Charge -- look at atoms that are only part of the bead (no ring rule)
    chg = 0
    for i in partitioning.keys():
        if partitioning[i] == cg_bead:
            chg += molecule.GetAtomWithIdx(i).GetFormalCharge()
    # Chem.rdmolfiles.SmilesWriter.write(sw,frag.GetMol())
    Chem.rdmolfiles.SmilesWriter.write(sw, Chem.rdmolops.AddHs(frag.GetMol(), addCoords=True))
    sw.close()
    # Read file
    s = ''
    try:
        f = open(tmpfile, 'r')
        s = f.readlines()
        f.close()
    except IOError as e:
        print("Error. Can't read file", tmpfile)
        print(e)
        exit(1)
    os.remove(tmpfile)
    return s[1].split()[0], wc_log_p, chg


def smi2alogps(smi, wc_log_p, bead, trial=False):
    """Returns water/octanol partitioning free energy
    according to ALOGPS"""
    session = requests.session()
    req = session.get('http://vcclab.org/web/alogps/calc?SMILES=' + str(smi))
    doc = BeautifulSoup.BeautifulSoup(req.content)
    soup = doc.prettify()
    found_mol_1 = False
    log_p = ''
    for line in soup.split("\n"):
        line = line.split()
        if "mol_1" in line:
            log_p = float(line[line.index('mol_1') + 1])
            found_mol_1 = True
            break
    if not found_mol_1:
        # If we're forcing a prediction, use Wildman-Crippen
        if args.forcepred:
            if trial:
                wrn = "; Warning: bead ID " + str(bead) + \
                      " predicted from Wildman-Crippen. Fragment " + str(smi) + "\n"
                sys.stderr.write(wrn)
            log_p = wc_log_p
        else:
            if args.verbose:
                print("ALOGPS can't predict fragment:", smi)
            exit(1)
    return convert_log_k(log_p)


def convert_log_k(log_k):
    """Convert log_{10}K to free energy (in kJ/mol)"""
    return 0.008314 * 300.0 * log_k / math.log10(math.exp(1))


def mad(bead_type, delta_f, in_ring=False):
    """Mean absolute difference between type type and delta_f"""
    if in_ring:
        # 3 beads in ring have the same atom type. Their
        # sum needs to reproduce the free energy.
        return math.fabs(3 * delta_f_types[bead_type] - 3 * delta_f)
    else:
        return math.fabs(delta_f_types[bead_type] - delta_f)


def determine_bead_type(delta_f, charge, hbonda, hbondd, in_ring):
    """Determine CG bead type from delta_f value, charge,
    and hbond acceptor, and donor"""
    # if args.verbose:
    # print "; dF:",delta_f/4.2,'kcal/mol'
    if charge < -1 or charge > +1:
        print("Charge is too large:", charge)
        print("No adequate force-field parameter")
        exit(1)
    bead_type = []
    if in_ring:
        # We're in a ring, divide free energy by 3
        # (average number of beads per ring)
        if abs(delta_f) > 0.1:
            delta_f *= 2 / 3.
    if charge != 0:
        # The compound has a +/- charge -> Q type
        if hbonda > 0 and hbondd > 0:
            bead_type = "Qda"
        elif hbonda > 0 and hbondd == 0:
            bead_type = "Qa"
        elif hbonda == 0 and hbondd > 0:
            bead_type = "Qd"
        else:
            bead_type = "Q0"
    else:
        # Neutral group
        # Use Hbond information only if we're close to Nda, Na, Nd types
        error = mad('Nda', delta_f, in_ring)
        if error < 2.0 and (hbonda > 0 or hbondd > 0):
            if hbonda > 0 and hbondd > 0:
                bead_type = "Nda"
            elif hbonda > 0 and hbondd == 0:
                bead_type = "Na"
            elif hbonda == 0 and hbondd > 0:
                bead_type = "Nd"
        else:
            # all other cases. Simply find the atom type that's closest in
            # free energy.
            othertypes = ['P5', 'P4', 'P3', 'P2', 'P1', 'N0', 'C5', 'C4', 'C3', 'C2', 'C1']
            min_error = 1000.0
            for cgtype in othertypes:
                tmp_error = mad(cgtype, delta_f, in_ring)
                if tmp_error < min_error:
                    min_error = tmp_error
                    bead_type = cgtype
    # if error > 5:
    # print "Warning: large error between beead type and log_k value:"
    # print " type {:s} ({:5.2f}) vs. {:5.2f}".format(bead_type,
    # delta_f_types[bead_type],delta_f)
    if in_ring:
        bead_type = "S" + bead_type
    return bead_type


def check_additivity(beadtypes, molecule):
    """Check additivity assumption between sum of free energies of CG beads
    and free energy of whole molecule"""
    # If there's only one bead, don't check.
    if len(beadtypes) == 1:
        return True
    sum_frag = 0.0
    rings = False
    for bead in beadtypes:
        if bead[0] == "S":
            bead = bead[1:]
            rings = True
        sum_frag += delta_f_types[bead]
    # Wildman-Crippen log_p
    wc_log_p = rdMolDescriptors.CalcCrippenDescriptors(molecule)[0]
    # Write out SMILES string of entire molecule
    tmpfile = 'tmp-auto-martini.smi'
    sw = rdmolfiles.SmilesWriter(tmpfile, nameHeader='')
    rdmolfiles.SmilesWriter.write(sw, molecule)
    sw.close()
    # Read file
    s = ''
    try:
        f = open(tmpfile, 'r')
        s = f.readlines()
        f.close()
    except IOError as e:
        print("Error. Can't read file", tmpfile)
        print(e)
        exit(1)
    os.remove(tmpfile)
    whole_mol_dg = smi2alogps(s[1].split()[0], wc_log_p, "MOL", True)
    m_ad = math.fabs((whole_mol_dg - sum_frag) / whole_mol_dg)
    if args.verbose:
        print("; Mapping additivity assumption ratio: {0:7.4f} ({1:7.4f} vs {2:7.4f})".format(m_ad,
                                                                                              whole_mol_dg, sum_frag))
    if (not rings and m_ad < 0.5) or rings:
        return True
    else:
        return False


def print_atoms(arguments, cgbeads, molecule, hbonda, hbondd, partitioning, ringatoms, ringatoms_flat, trial=False):
    """print CG Atoms in itp format"""
    atomnames = []
    beadtypes = []

    for bead in range(len(cgbeads)):
        # Determine SMI of substructure
        try:
            smi_frag, wc_log_p, charge = substruct2smi(molecule, partitioning, bead, cgbeads, ringatoms)
        except (NameError, TypeError, ValueError):
            return []
        atom_name = ""
        for character, count in sorted(six.iteritems(letter_occurrences(smi_frag))):
            try:
                float(character)
                pass
            except ValueError:
                if count == 1:
                    atom_name += "{:s}".format(character)
                else:
                    atom_name += "{:s}{:s}".format(character, str(count))
        # Get charge for smi_frag
        mol_frag = gen_molecule_smi(smi_frag)
        charge_frag = get_charge(mol_frag)
        # Extract ALOGPS free energy
        try:
            if charge_frag == 0:
                alogps = smi2alogps(smi_frag, wc_log_p, bead + 1, trial)
            else:
                alogps = 0.0
        except (NameError, TypeError, ValueError):
            return []
        hbond_a_flag = 0
        for at in hbonda:
            if partitioning[at] == bead and at != 0:
                hbond_a_flag = 1
                break
        hbond_d_flag = 0
        for at in hbondd:
            if partitioning[at] == bead and at != 0:
                hbond_d_flag = 1
                break
        in_ring = cgbeads[bead] in ringatoms_flat
        bead_type = determine_bead_type(alogps, charge, hbond_a_flag, hbond_d_flag, in_ring)
        atom_name = ""
        name_index = 0
        while atom_name in atomnames or name_index == 0:
            name_index += 1
            atom_name = "{:1s}{:02d}".format(bead_type[0], name_index)
        atomnames.append(atom_name)
        if not trial:
            print('    {:<5d} {:5s}     1             {:5s}     {:7s} {:<5d}    {:2d}         ; {:s}'.format(
                bead + 1, bead_type, arguments.molname, atom_name, bead + 1, charge, smi_frag))
        beadtypes.append(bead_type)
    return atomnames, beadtypes


def print_bonds(cgbeads, molecule, partitioning, cgbead_coords, ringatoms, trial=False):
    """print CG bonds in itp format"""
    if not trial:
        print("")
    # Bond information
    bondlist = []
    constlist = []
    if len(cgbeads) > 1:
        for i in range(len(cgbeads)):
            for j in range(i + 1, len(cgbeads)):
                dist = np.linalg.norm(cgbead_coords[i] - cgbead_coords[j]) * 0.1
                if dist < 0.65:
                    # Are atoms part of the same ring
                    in_ring = False
                    for ring in ringatoms:
                        if cgbeads[i] in ring and cgbeads[j] in ring:
                            in_ring = True
                            break
                    if in_ring:
                        constlist.append([i, j, dist])
                    else:
                        # Check that the bond is not too short
                        if dist < 0.2:
                            raise NameError('Bond too short')
                        # Look for a bond between an atom of i and an atom of j
                        found_connection = False
                        atoms_in_bead_i = []
                        for aa in partitioning.keys():
                            if partitioning[aa] == i:
                                atoms_in_bead_i.append(aa)
                        atoms_in_bead_j = []
                        for aa in partitioning.keys():
                            if partitioning[aa] == j:
                                atoms_in_bead_j.append(aa)
                        for ib in range(len(molecule.GetBonds())):
                            abond = molecule.GetBondWithIdx(ib)
                            if (abond.GetBeginAtomIdx() in atoms_in_bead_i and
                                abond.GetEndAtomIdx() in atoms_in_bead_j) or \
                                    (abond.GetBeginAtomIdx() in atoms_in_bead_j and
                                     abond.GetEndAtomIdx() in atoms_in_bead_i):
                                found_connection = True
                        if found_connection:
                            bondlist.append([i, j, dist])

        for ring in ringatoms:
            # Only keep one bond between a ring and a given external bead
            for i in range(len(cgbeads)):
                at = cgbeads[i]
                if at not in ring:
                    bonds_to_ring = []
                    for b in bondlist:
                        if (b[0] in ring and b[1] == at) or \
                                (b[0] == at and b[1] in ring):
                            bonds_to_ring.append(b)
                    # keep closest
                    closest_bond = [-1, -1, 1000.0]
                    for r in range(len(bonds_to_ring)):
                        if bonds_to_ring[r][2] < closest_bond[2]:
                            closest_bond = bonds_to_ring[r]
                    # Delete the other bonds
                    for b in bonds_to_ring:
                        if b != closest_bond:
                            bondlist.remove(b)
            bead_bonded_to_ring = []
            for i in range(len(cgbeads)):
                atoms_in_bead = []
                for key, val in six.iteritems(partitioning):
                    if val == i:
                        atoms_in_bead.append(key)
                for b in bondlist:
                    if (b[0] in ring and b[1] in atoms_in_bead) or \
                            (b[0] in atoms_in_bead and b[1] in ring):
                        bead_bonded_to_ring.append(i)
            # Delete bond between 2 beads if they're both linked
            # to the same ring.
            for i in range(len(cgbeads)):
                for j in range(i + 1, len(cgbeads)):
                    if cgbeads[i] in bead_bonded_to_ring and \
                       cgbeads[j] in bead_bonded_to_ring:
                        for b in bondlist:
                            if (b[0] == i and b[1] == j) or \
                                    (b[0] == j and b[1] == i):
                                bondlist.remove(b)

        # Replace bond by constraint if both atoms have constraints
        # to the same third atom
        bond_list_idx = 0
        while bond_list_idx < len(bondlist):
            b = bondlist[bond_list_idx]
            for i in range(len(cgbeads)):
                const_i = False
                const_j = False
                for c in constlist:
                    if (c[0] == b[0] and c[1] == i) or \
                            (c[0] == i and c[1] == b[0]):
                        const_i = True
                    if (c[0] == b[1] and c[1] == i) or \
                            (c[0] == i and c[1] == b[1]):
                        const_j = True
                if const_i and const_j:
                    constlist.append(b)
                    bondlist.remove(b)
            bond_list_idx += 1

        # Go through list of constraints. If we find an extra
        # possible constraint between beads that have constraints,
        # add it.
        beads_with_const = []
        for c in constlist:
            if c[0] not in beads_with_const:
                beads_with_const.append(c[0])
            if c[1] not in beads_with_const:
                beads_with_const.append(c[1])
        beads_with_const = sorted(beads_with_const)
        for i in range(len(beads_with_const)):
            for j in range(1 + i, len(beads_with_const)):
                const_exists = False
                for c in constlist:
                    if (c[0] == i and c[1] == j) or \
                            (c[0] == j and c[1] == i):
                        const_exists = True
                        break
                if not const_exists:
                    dist = np.linalg.norm(cgbead_coords[i] - cgbead_coords[j]) * 0.1
                    if dist < 0.35:
                        # Check that it's not in the bond list
                        in_bond_list = False
                        for b in bondlist:
                            if (b[0] == i and b[1] == j) or \
                                    (b[0] == j and b[0] == i):
                                in_bond_list = True
                                break
                        if not in_bond_list:
                            constlist.append([i, j, dist])

        if not trial:
            if len(bondlist) > 0:
                print("[bonds]")
                print("; i j     funct     length    force.c.")
                for b in bondlist:
                    # Make sure atoms in bond are not part of the same ring
                    print("    {:d} {:d}     1             {:4.2f}        1250".format(
                        b[0] + 1, b[1] + 1, b[2]))
                print("")
            if len(constlist) > 0:
                print("[constraints]")
                print(";  i   j     funct   length")
                for c in constlist:
                    print("   {:<3d} {:<3d}   1       {:4.2f}".format(
                        c[0] + 1, c[1] + 1, c[2]))
                print("")
            # Make sure there's at least a bond to every atom
            for i in range(len(cgbeads)):
                bond_to_i = False
                for b in bondlist + constlist:
                    if i in [b[0], b[1]]:
                        bond_to_i = True
                if not bond_to_i:
                    print("Error. No bond to atom", i + 1)
                    exit(1)
    return bondlist, constlist


def print_angles(cgbeads, molecule, partitioning, cgbead_coords, bondlist, constlist, ringatoms):
    """print CG angles in itp format"""
    if len(cgbeads) > 2:
        # Angles
        angle_list = []
        for i in range(len(cgbeads)):
            for j in range(len(cgbeads)):
                for k in range(len(cgbeads)):
                    # Only up to 2 atoms can be ring-like
                    all_in_ring = False
                    for ring in ringatoms:
                        if cgbeads[i] in ring and cgbeads[j] in ring and \
                           cgbeads[k] in ring:
                            all_in_ring = True
                            break
                    # Forbid all atoms linked by constraints
                    all_constraints = False
                    ij_bonded = False
                    jk_bonded = False
                    ij_const = False
                    jk_const = False
                    for b in bondlist + constlist:
                        if i in [b[0], b[1]] and j in [b[0], b[1]]:
                            ij_bonded = True
                            if b in constlist:
                                ij_const = True
                        if j in [b[0], b[1]] and k in [b[0], b[1]]:
                            jk_bonded = True
                            if b in constlist:
                                jk_const = True
                    if ij_const and jk_const:
                        all_constraints = True
                    if not all_in_ring and ij_bonded and jk_bonded and \
                       i != j and j != k and i != k and \
                            not all_constraints:
                        # Measure angle between i, j, and k.
                        angle = 180. / math.pi * math.acos(
                            np.dot(cgbead_coords[i] - cgbead_coords[j],
                                   cgbead_coords[k] - cgbead_coords[j]) / (
                                np.linalg.norm(cgbead_coords[i] - cgbead_coords[j]) *
                                np.linalg.norm(cgbead_coords[k] - cgbead_coords[j])))
                        # Look for any double bond between atoms belonging to these CG beads.
                        atoms_in_fragment = []
                        for aa in partitioning.keys():
                            if partitioning[aa] == j:
                                atoms_in_fragment.append(aa)
                        forc_const = 25.0
                        for ib in range(len(molecule.GetBonds())):
                            abond = molecule.GetBondWithIdx(ib)
                            if abond.GetBeginAtomIdx() in atoms_in_fragment and \
                               abond.GetEndAtomIdx() in atoms_in_fragment:
                                bondtype = molecule.GetBondBetweenAtoms(abond.GetBeginAtomIdx(),
                                                                        abond.GetEndAtomIdx()).GetBondType()
                                if bondtype == rdchem.BondType.DOUBLE:
                                    forc_const = 45.0
                        new_angle = True
                        for a in angle_list:
                            if i in a and j in a and k in a:
                                new_angle = False
                        if new_angle:
                            angle_list.append([i, j, k, angle, forc_const])
        if len(angle_list) > 0:
            print("[angles]")
            print("; i j k         funct   angle   force.c.")
            for a in angle_list:
                print("  {:d} {:d} {:d}         2       {:<5.1f}  {:5.1f}".format(
                    a[0] + 1, a[1] + 1, a[2] + 1, a[3], a[4]))
            print("")
    return


def print_dihedrals(cgbeads, constlist, ringatoms, cgbead_coords):
    """Print CG dihedrals in itp format"""
    if len(cgbeads) > 3:
        # Dihedrals
        dihed_list = []
        # Three ring atoms and one non ring
        for i in range(len(cgbeads)):
            for j in range(len(cgbeads)):
                for k in range(len(cgbeads)):
                    for l in range(len(cgbeads)):
                        if i != j and i != k and i != l \
                                and j != k and j != l and k != l:
                            # 3 atoms need to be ring like (in one ring!)
                            three_in_ring = False
                            for ring in ringatoms:
                                if [[cgbeads[i] in ring], [cgbeads[j] in ring],
                                   [cgbeads[k] in ring], [cgbeads[l] in ring]].count([True]) >= 3:
                                    three_in_ring = True
                                    break
                            for b in constlist:
                                if i in [b[0], b[1]] and j in [b[0], b[1]]:
                                    pass
                                if j in [b[0], b[1]] and k in [b[0], b[1]]:
                                    pass
                                if k in [b[0], b[1]] and l in [b[0], b[1]]:
                                    pass
                            # Distance criterion--beads can't be far apart
                            disthres = 0.35
                            close_enough = False
                            if np.linalg.norm(cgbead_coords[i] - cgbead_coords[j]) * 0.1 < disthres \
                                    and np.linalg.norm(cgbead_coords[j] - cgbead_coords[k]) * 0.1 < disthres \
                                    and np.linalg.norm(cgbead_coords[k] - cgbead_coords[l]) * 0.1 < disthres:
                                close_enough = True
                            already_dih = False
                            for dih in dihed_list:
                                if dih[0] == l and dih[1] == k and \
                                   dih[2] == j and dih[3] == i:
                                    already_dih = True
                                    break
                            if three_in_ring and close_enough and not already_dih:
                                r1 = cgbead_coords[j] - cgbead_coords[i]
                                r2 = cgbead_coords[k] - cgbead_coords[j]
                                r3 = cgbead_coords[l] - cgbead_coords[k]
                                p1 = np.cross(r1, r2) / (np.linalg.norm(r1) *
                                                         np.linalg.norm(r2))
                                p2 = np.cross(r2, r3) / (np.linalg.norm(r2) *
                                                         np.linalg.norm(r3))
                                r2 /= np.linalg.norm(r2)
                                cosphi = np.dot(p1, p2)
                                sinphi = np.dot(r2, np.cross(p1, p2))
                                angle = 180. / math.pi * arctan2(sinphi, cosphi)
                                forc_const = 10.0
                                dihed_list.append([i, j, k, l, angle, forc_const])
        if len(dihed_list) > 0:
            print("[dihedrals]")
            print(";  i     j    k    l   funct   angle  force.c.")
            for d in dihed_list:
                print("   {:d}     {:d}    {:d}    {:d}       2     {:<5.1f}  {:5.1f}".format(
                    d[0] + 1, d[1] + 1, d[2] + 1, d[3] + 1, d[4], d[5]))
            print("")
    return


if __name__ == '__main__':

    # Parse command-line options
    parser = argparse.ArgumentParser(description='Map atomistic structure to MARTINI mapping',
                                     epilog='Tristan BEREAU (2014)')
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--sdf', dest='sdf', type=str, required=False, help='SDF file of atomistic coordinates')
    group.add_argument('--smi', dest='smi', type=str, required=False, help='SMILES string of atomistic structure')
    parser.add_argument('--mol', dest='molname', type=str, required=True, help='Name of CG molecule')
    parser.add_argument('--xyz', dest='xyz', type=str, help='output atomistic structure to OUT.xyz file')
    parser.add_argument('--gro', dest='gro', type=str, help='output CG structure to OUT.gro file')
    parser.add_argument('--verbose', dest='verbose', action='store_true', help='verbose')
    parser.add_argument('--fpred', dest='forcepred', action='store_true', help='verbose')

    args = parser.parse_args()

    if args.sdf:
        # Generate molecule's structure from SDF
        mol = gen_molecule_sdf(args.sdf)
    else:
        mol = gen_molecule_smi(args.smi)

    # Get molecule's features
    feats = extract_features(mol)

    # Identify ring-type atoms
    ring_atoms = get_ring_atoms(feats)

    # Get Hbond information
    hbond_a = get_hbond_a(feats)
    hbond_d = get_hbond_d(feats)

    # Flatten list of ring atoms
    ring_atoms_flat = list(chain.from_iterable(ring_atoms))

    # Optimize CG bead positions -- keep all possibilities in case something goes
    # wrong later in the code.
    listCGBeads, listBeadPos = find_bead_pos(mol, ring_atoms)

    # Loop through best 1% cg_beads and avg_pos
    atom_names = []
    cg_bead_coords = []
    maxAttempts = int(math.ceil(0.5 * len(listCGBeads)))
    if args.verbose:
        print("; Max. number of attempts:", maxAttempts)
    attempt = 0
    while attempt < maxAttempts:
        cg_beads = listCGBeads[attempt]
        bead_pos = listBeadPos[attempt]
        success = True
        # Extract atom coordinates of heavy atoms
        heavy_atom_coords = get_heavy_atom_coords(mol)
        # Extract position of CG beads
        cg_bead_coords = get_cg_bead_coords(mol, cg_beads, bead_pos, ring_atoms_flat)

        # Partition atoms into CG beads
        atom_partitioning = voronoi_atoms(cg_bead_coords, heavy_atom_coords)
        if args.verbose:
            print("; Atom partitioning:", atom_partitioning)

        atom_names, bead_types = print_atoms(args, cg_beads, mol, hbond_a, hbond_d, atom_partitioning,
                                             ring_atoms, ring_atoms_flat, True)
        if not atom_names:
            success = False
        # Check additivity between fragments and entire molecule
        if not check_additivity(bead_types, mol):
            success = False
        # Bond list
        try:
            bond_list, const_list = print_bonds(cg_beads, mol, atom_partitioning, cg_bead_coords, ring_atoms, True)
        except (NameError, ValueError):
            success = False

        if success:
            print_header(args)
            atom_names, bead_types = print_atoms(args, cg_beads, mol, hbond_a, hbond_d, atom_partitioning,
                                                 ring_atoms, ring_atoms_flat, False)
            bond_list, const_list = print_bonds(cg_beads, mol, atom_partitioning, cg_bead_coords, ring_atoms, False)
            print_angles(cg_beads, mol, atom_partitioning, cg_bead_coords, bond_list, const_list, ring_atoms)
            print_dihedrals(cg_beads, const_list, ring_atoms, cg_bead_coords)
            # We've reached all the way here, exit the while loop
            attempt = maxAttempts + 1
        else:
            attempt += 1
    if attempt == maxAttempts:
        err = "; ERROR: no successful mapping found.\n" + \
              "; Try running with the '--fpred' and/or '--verbose' options.\n"
        sys.stderr.write(err)
        exit(1)

    # Optional atomistic output to XYZ file
    if args.xyz:
        output_xyz(args, mol)
    # Optional CG output to GRO file
    if args.gro:
        output_gro(args, cg_bead_coords, atom_names)
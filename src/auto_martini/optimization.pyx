# !/usr/bin/python
#  -*- coding: utf8 -*-

'''
  Created on March 13, 2019 by Andrew Abi-Mansour

          _    _ _______ ____    __  __          _____ _______ _____ _   _ _____ 
     /\  | |  | |__   __/ __ \  |  \/  |   /\   |  __ \__   __|_   _| \ | |_   _|
    /  \ | |  | |  | | | |  | | | \  / |  /  \  | |__) | | |    | | |  \| | | |  
   / /\ \| |  | |  | | | |  | | | |\/| | / /\ \ |  _  /  | |    | | | . ` | | |  
  / ____ \ |__| |  | | | |__| | | |  | |/ ____ \| | \ \  | |   _| |_| |\  |_| |_ 
 /_/    \_\____/   |_|  \____/  |_|  |_/_/    \_\_|  \_\ |_|  |_____|_| \_|_____|
                                                                                 
                                                                                 
  @originally written by: Tristan BEREAU (bereau at mpip-mainz.mpg.de)
  @modified by: Andrew Abi-Mansour (andrew.gaam at gmail.com)

  Auto_martini: a tool for automatic MARTINI mapping and parametrization of small organic molecules                                               
    
  Auto_martini is open-source, distributed under the terms of the GNU Public
  License, version 2 or later. It is distributed in the hope that it will
  be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
  of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. You should have
  received a copy of the GNU General Public License along with PyGran.
  If not, see http://www.gnu.org/licenses . See also top-level README
  and LICENSE files.

  Github link to original repo: https://github.com/tbereau/auto_martini
'''

from auto_martini.engine.common import *

import numpy as np

try:
  cimport numpy as np
except:
  pass # no Cython .. no problem

def read_bead_params():
    """Returns bead parameter dictionary
    CG Bead vdw radius (in Angstroem)"""
    bead_params = dict()
    bead_params['rvdw'] = 4.7 / 2.
    bead_params['rvdw_aromatic'] = 4.3 / 2.
    bead_params['rvdw_cross'] = 0.5*((4.7 / 2.) + (4.3 / 2.))
    bead_params['offset_bd_weight'] = 50.
    bead_params['offset_bd_aromatic_weight'] = 20.
    bead_params['lonely_atom_penalize'] = 0.28 #0.20
    bead_params['bd_bd_overlap_coeff'] = 9.0
    bead_params['at_in_bd_coeff'] = 0.9
    return bead_params

def gaussian_overlap(conformer, bead1, bead2, ringatoms):
    """"Returns overlap coefficient between two gaussians
    given distance dist"""
    logging.debug('Entering gaussian_overlap()')
    dist = Chem.rdMolTransforms.GetBondLength(conformer, int(bead1), int(bead2))
    bead_params = read_bead_params()
    sigma = bead_params['rvdw']
    if bead1 in ringatoms and bead2 in ringatoms:
        sigma = bead_params['rvdw_aromatic']
    if bead1 in ringatoms and bead2 not in ringatoms or \
       bead1 not in ringatoms and bead2 in ringatoms:
        sigma = bead_params['rvdw_cross']
    return bead_params['bd_bd_overlap_coeff'] * math.exp(-dist ** 2 / 4. / sigma ** 2)


def atoms_in_gaussian(molecule, conformer, bead_id, ringatoms):
    """Returns weighted sum of atoms contained in bead bead_id"""
    logging.debug('Entering atoms_in_gaussian()')
    weight_sum = 0.0
    bead_params = read_bead_params()
    sigma = bead_params['rvdw']
    lumped_atoms = []
    if bead_id in ringatoms:
        sigma = bead_params['rvdw_aromatic']
    for i in range(conformer.GetNumAtoms()):
        dist_bd_at = Chem.rdMolTransforms.GetBondLength(conformer, i, int(bead_id))
        if dist_bd_at < sigma:
            lumped_atoms.append(i)
        weight_sum -= molecule.GetAtomWithIdx(i).GetMass() * math.exp(-dist_bd_at ** 2 / 2 / sigma ** 2)
    return bead_params['at_in_bd_coeff'] * weight_sum, lumped_atoms


def penalize_lonely_atoms(molecule, conformer, lumped_atoms):
    """Penalizes configuration if atoms aren't included
    in any CG bead"""
    logging.debug('Entering penalize_lonely_atoms()')
    weight_sum = 0.0
    bead_params = read_bead_params()
    num_atoms = conformer.GetNumAtoms()
    atoms_array = np.arange(num_atoms)
    for i in np.nditer(np.arange(atoms_array.size)):
        if atoms_array[i] not in lumped_atoms:
            weight_sum += molecule.GetAtomWithIdx(int(atoms_array[i])).GetMass()
    return bead_params['lonely_atom_penalize'] * weight_sum


def eval_gaussian_interac(molecule, conformer, list_beads, ringatoms):
    """From collection of CG beads placed on mol, evaluate
    objective function of interacting beads"""
    logging.debug('Entering eval_gaussian_interac()')

    weight_sum = 0.0
    weight_overlap = 0.0
    weight_at_in_bd = 0.0
    bead_params = read_bead_params()
    
    # Offset energy for every new CG bead.
    # Distinguish between aromatics and others.
    num_aromatics = 0
    lumped_atoms = []

    # Creat list_beads array and loop over indeces
    list_beads_array = np.asarray(list_beads)
    for i in np.nditer(np.arange(list_beads_array.size)):
        if list_beads_array[i] in ringatoms:
            num_aromatics += 1
    weight_offset_bd_weights = bead_params['offset_bd_weight'] * (list_beads_array.size - num_aromatics) + \
        bead_params['offset_bd_aromatic_weight'] * num_aromatics
    weight_sum += weight_offset_bd_weights

    # Repulsive overlap between CG beads
    for i in np.nditer(np.arange(list_beads_array.size)):
        if i < list_beads_array.size-1:
            for j in np.nditer(np.arange(i+1, list_beads_array.size)):
                weight_overlap += gaussian_overlap(conformer, list_beads_array[i], list_beads_array[j], ringatoms)
    weight_sum += weight_overlap

    # Attraction between atoms nearby to CG bead
    for i in np.nditer(np.arange(list_beads_array.size)):
        weight, lumped = atoms_in_gaussian(molecule, conformer, list_beads_array[i], ringatoms)
        weight_at_in_bd += weight
        lumped_array = np.asarray(lumped)
        for j in np.nditer(np.arange(lumped_array.size)):
            if lumped_array[j] not in lumped_atoms:
                lumped_atoms.append(lumped_array[j])
    weight_sum += weight_at_in_bd
    # Penalty for excluding atoms
    weight_lonely_atoms = penalize_lonely_atoms(molecule, conformer, lumped_atoms)
    weight_sum += weight_lonely_atoms
    logging.debug(weight_sum, weight_offset_bd_weights, weight_overlap, weight_at_in_bd, weight_lonely_atoms)
    return weight_sum


def check_beads(list_heavyatoms, heavyatom_coords, trial_comb, ring_atoms, listbonds):
    """Check if CG bead positions in trailComb are acceptable"""
    logging.debug('Entering check_beads()')
    acceptable_trial = ''
    # Check for beads at the same place
    count = Counter(trial_comb)
    all_different = True
    for val in count.values():
        if val != 1:
            all_different = False
            acceptable_trial = False
            logging.debug('Error. Multiple beads on the same atom position for %s' % trial_comb)
            break
    if all_different:
        acceptable_trial = True
        # Check for beads linked by chemical bond (except in rings)
        bonds_in_rings = [0] * len(ring_atoms)
        for bi in range(len(trial_comb)):
            for bj in range(bi + 1, len(trial_comb)):
                if [trial_comb[bi], trial_comb[bj]] in listbonds \
                        or [trial_comb[bj], trial_comb[bi]] in listbonds:
                    bond_in_ring = False
                    for r in range(len(ring_atoms)):
                        if trial_comb[bi] in ring_atoms[r] and trial_comb[bj] in ring_atoms[r]:
                            bonds_in_rings[r] += 1
                            bond_in_ring = True
                    if not bond_in_ring:
                        acceptable_trial = False
                        logging.debug('Error. No bond in ring for %s' % trial_comb)
                        break
        if acceptable_trial:
            # Don't allow bonds between atoms of the same ring.
            for bir in range(len(bonds_in_rings)):
                if bonds_in_rings[bir] > 0:
                    logging.debug('Error. Bonds between atoms of the same ring for %s', trial_comb)
                    acceptable_trial = False
        if acceptable_trial:
            # Check for two terminal beads linked by only one atom
            for bi in range(len(trial_comb)):
                for bj in range(bi + 1, len(trial_comb)):
                    if ([item for sublist in listbonds for item in
                         sublist].count(trial_comb[bi]) == 1) and ([item for sublist
                                                                    in listbonds for item in sublist].count(
                            trial_comb[bj]) == 1):
                        # Both beads are on terminal atoms. Block contribution
                        # if the two terminal atoms are linked to the same atom.
                        partneri = ''
                        partnerj = ''
                        for bond in listbonds:
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
                            logging.debug('Error. Two terminal beads linked to the same atom for %s' % trial_comb)
    return acceptable_trial


def find_bead_pos(molecule, conformer, list_heavy_atoms, heavyatom_coords, ring_atoms, ringatoms_flat):
    """Try out all possible combinations of CG beads up to threshold number of beads per atom. Find
    arrangement with best energy score. Return all possible arrangements sorted by energy score."""
    
    logging.debug('Entering find_bead_pos()')
    
    # Check number of heavy atoms
    if len(list_heavy_atoms) == 0:
        print('Error. No heavy atom found.')
        exit(1)
    
    if len(list_heavy_atoms) == 1:
        # Put one CG bead on the one heavy atom.
        best_trial_comb = np.array(list(itertools.combinations(range(len(list_heavy_atoms)), 1)))
        avg_pos = [[conformer.GetAtomPosition(best_trial_comb[0])[j] for j in range(3)]]
        return best_trial_comb, avg_pos
    
    if len(list_heavy_atoms) > 50:
        print('Error. Exhaustive enumeration can\'t handle large molecules.')
        print('Number of heavy atoms: %d' % len(list_heavy_atoms))
        exit(1)
    # List of bonds between heavy atoms
    list_bonds = []
    for i in range(len(list_heavy_atoms)):
        for j in range(i + 1, len(list_heavy_atoms)):
            if molecule.GetBondBetweenAtoms(int(list_heavy_atoms[i]), int(list_heavy_atoms[j])) is not None:
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
    
    for num_beads in range(1, max_beads + 1):

        # Use recursive function to loop through all possible
        # combinations of CG bead positions.
        seq_one_beads = np.array(list(itertools.combinations(list_heavy_atoms, num_beads)))
        combs = []
        energies = []

        # Trial positions: any heavy atom
        for seq in seq_one_beads:
            trial_comb = list(seq)
            acceptable_trial = check_beads(list_heavy_atoms, heavyatom_coords, trial_comb, ring_atoms, list_bonds)
            if acceptable_trial:

                # Do the energy evaluation
                trial_ene = eval_gaussian_interac(molecule, conformer, trial_comb, ringatoms_flat)
                combs.append(trial_comb)
                energies.append(trial_ene)
                
                logging.info('; %s %s', trial_comb, trial_ene)
                # Make sure all atoms within one bead would be connected
                if all_atoms_in_beads_connected(trial_comb, heavyatom_coords, list_heavy_atoms, list_bonds):
                    # Accept the move
                    if trial_ene < ene_best_trial:
                        ene_best_trial = trial_ene
                        best_trial_comb = sorted(trial_comb)
                    # Get bead positions
                    beadpos = [[0]*3 for l in range(len(trial_comb))]
                    for l in range(len(trial_comb)):
                        beadpos[l] = [conformer.GetAtomPosition(int(sorted(trial_comb)[l]))[m] for m in range(3)]
                    # Store configuration
                    list_trial_comb.append([trial_comb, beadpos, trial_ene])

        if last_best_trial_comb == best_trial_comb:
            break

        last_best_trial_comb = best_trial_comb
        list_combs.append(combs)
        list_energies.append(energies)

    sorted_combs = np.array(sorted(list_trial_comb, key=itemgetter(2)))
    return sorted_combs[:, 0], sorted_combs[:, 1]

def all_atoms_in_beads_connected(trial_comb, heavyatom_coords, list_heavyatoms, bondlist):
    """Make sure all atoms within one CG bead are connected to at least
    one other atom in that bead"""
    logging.debug('Entering all_atoms_in_beads_connected()')
    # Bead coordinates are given by heavy atoms themselves
    cgbead_coords = []

    for i in range(len(trial_comb)):
        cgbead_coords.append(heavyatom_coords[list_heavyatoms.index(trial_comb[i])])
    voronoi = voronoi_atoms(cgbead_coords, heavyatom_coords)
    logging.debug('voronoi %s' % voronoi)

    for i in range(len(trial_comb)):
        cg_bead = trial_comb[i]
        num_atoms = list(voronoi.values()).count(voronoi[list_heavyatoms.index(cg_bead)])
        # sub-part of bond list that only contains atoms within CG bead
        sub_bond_list = []
        for j in range(len(bondlist)):
            if voronoi[list_heavyatoms.index(bondlist[j][0])] == voronoi[list_heavyatoms.index(cg_bead)] and \
               voronoi[list_heavyatoms.index(bondlist[j][1])] == voronoi[list_heavyatoms.index(cg_bead)]:
                sub_bond_list.append(bondlist[j])
        num_bonds = len(sub_bond_list)
        if num_bonds < num_atoms - 1 or num_atoms == 1:
            logging.debug('Error: Not all atoms in beads connected in %s' % trial_comb)
            logging.debug('Error: %s < %s, %s' % (num_bonds, num_atoms-1, sub_bond_list))
            return False
    return True


def voronoi_atoms(cgbead_coords, heavyatom_coords):
    """Partition all atoms between CG beads"""
    logging.debug('Entering voronoi_atoms()')
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
                logging.warning('Error. Can\'t find closest atom to bead %s' % i)
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
                    logging.warning('Error. Can\'t find an atom close to atom $s' % lonely_bead)
                    exit(1)
                partitioning[closest_bead] = lonely_bead
    return partitioning

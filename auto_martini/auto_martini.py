# !/usr/bin/python
#  -*- coding: utf8 -*-

'''
  Created on March 13, 2019

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

  BUGS found and fixed:

  - function substruct2smi() returned SMILES string in lower case letters
  - function printAtoms() uses undefined variables "hbondA" and "hbondD"
  - function checkAdditivity() defined "mad" variable which is already defined as a function
  - function genMoleculeSDF() did not use "sdf" input argument

'''

import sys,os
import argparse
import math
import numpy as np
from collections import Counter
from rdkit import Chem
from rdkit.Chem import ChemicalFeatures
from rdkit.Chem import rdchem
from rdkit.Chem import rdMolDescriptors
from rdkit import RDConfig
from itertools import chain
from operator import itemgetter
import lxml

import output
import topology

# For feature extraction
fdefName = os.path.join(RDConfig.RDDataDir,'BaseFeatures.fdef')
factory = ChemicalFeatures.BuildFeatureFactory(fdefName)

# Parameters
# CG Bead vdw radius (in Angstroem)
rvdw = 4.7/2.
rvdwAromatic = 4.3/2.
rvdwCross = 0.5*(rvdw+rvdwAromatic)
# Optimized parameters
offsetBeadWeight = 50
offsetBeadAromaticWeight = 20.
lonelyAtomPenalize = 0.20
bdBdOverlapCoeff = 9.0
atInBdCoeff = 0.9

def getHbondA(feats):
  '''Get Hbond acceptor information'''
  hbondA = []
  for feat in feats:
    if feat.GetFamily() == "Acceptor":
      for i in feat.GetAtomIds():
        if i not in hbondA:
          hbondA.append(i)
  return hbondA

def getHbondD(mol):
  '''Get Hbond donor information'''
  hbondD = []
  for feat in feats:
    if feat.GetFamily() == "Donor":
      for i in feat.GetAtomIds():
        if i not in hbondD:
          hbondD.append(i)
  return hbondD

def extractFeatures(mol):
  '''Extract features of mol'''
  feats = factory.GetFeaturesForMol(mol)
  # if len(feats) == 0:
  #   print "Error. Can't extract molecular features."
  #   exit(1)
  return feats

def getRingAtoms(feats, args):
  '''Get ring atoms'''
  ringAtoms = []
  for feat in feats:
    if feat.GetType() in ["RH6_6","RH5_5","RH4_4","RH3_3",\
      "Arom5","Arom6","Arom7","Arom8"]:
      newRing = []
      for at in feat.GetAtomIds():
          newRing.append(at)
      if newRing not in ringAtoms:
        ringAtoms.append(newRing)
  if args.verbose:
    print("; ring atoms:",ringAtoms)
  return ringAtoms

def gaussianOverlap(bead1, bead2, ringAtoms):
  '''Returns overlap coefficient between two gaussians
  given distance dist'''
  conf = mol.GetConformer()
  dist = Chem.rdMolTransforms.GetBondLength(conf,bead1,bead2)
  sigma = rvdw
  if bead1 in ringAtoms and bead2 in ringAtoms:
    sigma = rvdwAromatic
  if bead1 in ringAtoms and bead2 not in ringAtoms or \
    bead1 not in ringAtoms and bead2 in ringAtoms:
    sigma = rvdwCross
  return bdBdOverlapCoeff*math.exp(-dist**2/4./sigma**2)

def atomsInGaussian(mol, beadId, ringAtoms):
  '''Returns weighted sum of atoms contained in bead beadId'''
  weightSum = 0.0
  conf = mol.GetConformer()
  sigma = rvdw
  lumpedAtoms = []
  if beadId in ringAtoms:
    sigma = rvdwAromatic
  for i in range(conf.GetNumAtoms()):
    distBdAt = Chem.rdMolTransforms.GetBondLength(conf,i,beadId)
    if distBdAt < sigma:
      lumpedAtoms.append(i)
    weightSum -= mol.GetAtomWithIdx(i).GetMass() * \
      math.exp(-distBdAt**2/2/sigma**2)
  return atInBdCoeff*weightSum,lumpedAtoms

def penalizeLonelyAtoms(mol, lumpedAtoms):
  '''Penalizes configuration if atoms aren't included
  in any CG bead'''
  weightSum = 0.0
  conf = mol.GetConformer()
  for i in range(conf.GetNumAtoms()):
    if i not in lumpedAtoms:
      weightSum += mol.GetAtomWithIdx(i).GetMass()
  return lonelyAtomPenalize*weightSum

def evalGaussianInterac(mol, listBeads, ringAtoms):
  '''From collection of CG beads placed on mol, evaluate
  objective function of interacting beads'''
  weightSum = 0.0
  # Offset energy for every new CG bead.
  # Distinguish between aromatics and others.
  numAromatics = 0
  lumpedAtoms = []
  for i in listBeads:
    if i in ringAtoms:
      numAromatics += 1
  weightSum += offsetBeadWeight * (len(listBeads)-numAromatics) + \
    offsetBeadAromaticWeight * numAromatics
  # Repulsive overlap between CG beads
  for i in range(len(listBeads)):
    for j in range(i+1,len(listBeads)):
      weightSum += gaussianOverlap(listBeads[i], listBeads[j], ringAtoms)
  # Attraction between atoms nearby to CG bead
  for i in range(len(listBeads)):
    weight,lumped = atomsInGaussian(mol, \
      listBeads[i], ringAtoms)
    weightSum += weight
    for j in lumped:
      if j not in lumpedAtoms:
        lumpedAtoms.append(j)
  # Penalty for excluding atoms
  weightSum += penalizeLonelyAtoms(mol, lumpedAtoms)
  return weightSum

def enumerateSeq(heavyAtoms,depth):
  '''Enumerate all sequences of length depth among heavy atoms'''
  seq = [[heavyAtoms[0]] * depth]
  count = 0
  while seq[-1] != [heavyAtoms[-1]] * depth:
    lastEle = seq[-1]
    newEle = []
    lastIndex = 1
    toggleNext = True
    while lastIndex <= depth:
      if toggleNext:
        if lastEle[-lastIndex] == heavyAtoms[-1]:
          # find left-most number that's not == heavyAtoms[-1]
          lix = 0
          for lix in range(1,depth+1):
            if lastEle[-lix] != heavyAtoms[-1]:
              break
          for lx in range(0,lix-1):
            newEle = [heavyAtoms[heavyAtoms.index(lastEle[-lix]+1)]] + newEle
          lastIndex = lix - 1
          toggleNext = True
        else:
          newEle = [heavyAtoms[heavyAtoms.index(
            lastEle[-lastIndex])+1]] + newEle
          toggleNext = False
      else:
        newEle = [lastEle[-lastIndex]] + newEle
      lastIndex += 1
    if not toggleNext:
      seq.append(newEle)
  return seq

def findBeadPos(mol, feats, ringAtoms, args):
  '''Try out all possible combinations of CG beads
  up to threshold number of beads per atom. find
  arrangement with best energy score. Return all
  possible arrangements sorted by energy score.'''
  # number of atoms in mol
  conf = mol.GetConformer()
  numAtoms = conf.GetNumAtoms()
  # List of heavy atoms
  listHeavyAtoms = []
  if args.verbose:
    print("; Heavy atoms:"),
  for i in range(numAtoms):
    if mol.GetAtomWithIdx(i).GetSymbol() != "H":
      listHeavyAtoms.append(i)
      if args.verbose:
        print(mol.GetAtomWithIdx(i).GetSymbol()),
  if args.verbose:
    print("")
  if len(listHeavyAtoms) == 0:
    print("Error. No heavy atom found.")
    exit(1)
  if len(listHeavyAtoms) == 1:
    # Put one CG bead on the one heavy atom.
    bestTrialComb = enumerateSeq(listHeavyAtoms,1)[0]
    avgPos = [[conf.GetAtomPosition(bestTrialComb[0])[j] for j in range(3)]]
    return bestTrialComb,avgPos
  if len(listHeavyAtoms) > 25:
    print("Error. Exhaustive enumeration can't handle large molecules.")
    print("Number of heavy atoms:",len(listHeavyAtoms))
    exit(1)
  ringAtomsFlat = list(chain.from_iterable(ringAtoms))
  # List of bonds between heavy atoms
  listBonds = []
  for i in range(len(listHeavyAtoms)):
    for j in range(i+1,len(listHeavyAtoms)):
      if mol.GetBondBetweenAtoms(listHeavyAtoms[i],listHeavyAtoms[j]) != None:
        listBonds.append([listHeavyAtoms[i],listHeavyAtoms[j]])
  # Max number of beads. At most 2.5 heavy atoms per bead.
  maxBeads = int(len(listHeavyAtoms)/2.)
  # Collect all possible combinations of bead positions
  bestTrialComb = []
  listTrialComb = []
  eneBestTrial  = 1e6
  lastBestTrialComb = []

  # Keep track of all combinations and scores
  listCombs = []
  listEnergies = []

  # Heavy atom coordinates
  heavyAtomCoords = getHeavyAtomCoords(mol)

  for numBeads in range(1,maxBeads+1):
    # Use recursive function to loop through all possible
    # combinations of CG bead positions.
    seqOneBeads = enumerateSeq(listHeavyAtoms,numBeads)

    combs = []
    energies = []

    # Trial positions: any heavy atom
    for i in range(len(seqOneBeads)):
      trialComb = seqOneBeads[i]
      # Check for beads at the same place
      count = Counter(trialComb)
      allDifferent = True
      for val in count.values():
        if val != 1:
          allDifferent = False
          break
      if allDifferent:
        acceptableTrial = True
        # Check for beads linked by chemical bond (except in rings)
        acceptableTrial = True
        bondsInRings = [0]*len(ringAtoms)
        for bi in range(len(trialComb)):
          for bj in range(bi+1,len(trialComb)):
            if [trialComb[bi],trialComb[bj]] in listBonds \
              or [trialComb[bj],trialComb[bi]] in listBonds:
              bondInRing = False
              for r in range(len(ringAtoms)):
                if trialComb[bi] in ringAtoms[r] and trialComb[bj] in ringAtoms[r]:
                  bondsInRings[r] += 1
                  bondInRing = True
              if not bondInRing:
                acceptableTrial = False
                break
        if acceptableTrial:
          # Don't allow bonds between atoms of the same ring.
          for bir in range(len(bondsInRings)):
            if bondsInRings[bir] > 0:
              acceptableTrial = False
        if acceptableTrial:
          # Check for two terminal beads linked by only one atom
          for bi in range(len(trialComb)):
            for bj in range(bi+1,len(trialComb)):
              if ([item for sublist in listBonds for item in \
                sublist].count(trialComb[bi]) == 1) and ([item for sublist \
                in listBonds for item in sublist].count(trialComb[bj]) == 1):
                # Both beads are on terminal atoms. Block contribution
                # if the two terminal atoms are linked to the same atom.
                partneri = ''
                partnerj = ''
                for bond in listBonds:
                  if bond[0] == trialComb[bi]:
                    partneri = bond[1]
                  if bond[1] == trialComb[bi]:
                    partneri = bond[0]
                  if bond[0] == trialComb[bj]:
                    partnerj = bond[1]
                  if bond[1] == trialComb[bj]:
                    partnerj = bond[0]
                if partneri == partnerj:
                  acceptableTrial = False
                  break
        if acceptableTrial:
          # Don't accept single atom off a ring
          for bi in range(len(trialComb)):
            for r in ringAtomsFlat:
              bond = mol.GetBondBetweenAtoms(trialComb[bi], r)
              if bond is not None and \
                len(mol.GetAtomWithIdx(trialComb[bi]).GetBonds()) == 1:
                acceptableTrial = False
                break
        if acceptableTrial:
          # Do the energy evaluation
          trialEne  = evalGaussianInterac(mol, trialComb, ringAtomsFlat)
          combs.append(trialComb)
          energies.append(trialEne)
          if args.verbose:
            print(";",trialComb,trialEne)
          # Make sure all atoms within one bead would be connected
          if allAtomsInBeadsConnected(trialComb,
            heavyAtomCoords, listHeavyAtoms, listBonds):
            # Accept the move
            if trialEne < eneBestTrial:
              eneBestTrial  = trialEne
              bestTrialComb = sorted(trialComb)
            # Get bead positions
            beadPos = [[0]*3 for l in range(len(trialComb))]
            for l in range(len(trialComb)):
              beadPos[l] = [conf.GetAtomPosition(sorted(trialComb)[l])[m] for m in range(3)]
            # Store configuration
            listTrialComb.append([trialComb,beadPos,trialEne])
    if lastBestTrialComb == bestTrialComb:
      break
    lastBestTrialComb = bestTrialComb
    listCombs.append(combs)
    listEnergies.append(energies)
  if args.verbose:
    for at in bestTrialComb:
      print("; CG bead:",at)
    print(";      with energy:", eneBestTrial)
  sortedCombs = np.array(sorted(listTrialComb, key=itemgetter(2)))
  return sortedCombs[:,0],sortedCombs[:,1]

def getHeavyAtomCoords(mol):
  '''Extract atomic coordinates of heavy atoms in molecule mol'''
  heavyAtomCoords = []
  conf = mol.GetConformer()
  # number of atoms in mol
  numAtoms = mol.GetConformer().GetNumAtoms()
  for i in range(numAtoms):
    if mol.GetAtomWithIdx(i).GetSymbol() != "H":
      heavyAtomCoords.append(np.array(
        [conf.GetAtomPosition(i)[j] for j in range(3)]))
  return heavyAtomCoords

def getCGBeadCoords(mol, cgBeads, avgPos, ringAtomsFlat):
  '''Extract coordinates of CG beads'''
  # CG beads are averaged over best trial combinations for all
  # non-aromatic atoms.
  cgBeadCoords = []
  conf = mol.GetConformer()
  for i in range(len(cgBeads)):
    if cgBeads[i] in ringAtomsFlat:
      cgBeadCoords.append(np.array([conf.GetAtomPosition(cgBeads[i])[j]
        for j in range(3)]))
    else:
      # Use average
      cgBeadCoords.append(np.array(avgPos[i]))
  return cgBeadCoords

def allAtomsInBeadsConnected(trialComb, heavyAtomCoords, listHeavyAtoms, bondList):
  '''Make sure all atoms within one CG bead are connected to at least
  one other atom in that bead'''
  # Bead coordinates are given by heavy atoms themselves
  cgBeadCoords = []
  for i in range(len(trialComb)):
    cgBeadCoords.append(heavyAtomCoords[listHeavyAtoms.index(trialComb[i])])
  voronoi = voronoiAtoms(cgBeadCoords,heavyAtomCoords)
  for i in range(len(trialComb)):
    cgBead = trialComb[i]
    numAtoms = list(voronoi.values()).count(voronoi[listHeavyAtoms.index(cgBead)])
    # sub-part of bond list that only contains atoms within CG bead
    subBondList = []
    for j in range(len(bondList)):
      if voronoi[listHeavyAtoms.index(bondList[j][0])] == voronoi[listHeavyAtoms.index(cgBead)] and \
        voronoi[listHeavyAtoms.index(bondList[j][1])] == voronoi[listHeavyAtoms.index(cgBead)]:
        subBondList.append(bondList[j])
    numBonds = len(subBondList)
    if numBonds < numAtoms-1 or numAtoms == 1:
      return False
  return True

def voronoiAtoms(cgBeadCoords, heavyAtomCoords):
  '''Partition all atoms between CG beads'''
  atomPartitioning = {}
  for j in range(len(heavyAtomCoords)):
    if j not in atomPartitioning.keys():
      # Voronoi to check whether atom is closest to bead
      closestBead = True
      beadAt     = -1
      distBeadAt = 1000
      for k in range(len(cgBeadCoords)):
        distk = np.linalg.norm(cgBeadCoords[k]-heavyAtomCoords[j])
        if distk < distBeadAt:
          distBeadAt = distk
          beadAt = k
      atomPartitioning[j] = beadAt
  if len(cgBeadCoords) > 1:
    # Book-keeping of closest atoms to every bead
    closestAtoms = {}
    for i in range(len(cgBeadCoords)):
      closestAtom = -1
      closestDist = 10000.0
      for j in range(len(heavyAtomCoords)):
        distBeadAt = np.linalg.norm(cgBeadCoords[i]-heavyAtomCoords[j])
        if distBeadAt < closestDist:
          closestDist = distBeadAt
          closestAtom = j
      if closestAtom == -1:
        print("Error. Can't find closest atom to bead", i)
        exit(1)
      closestAtoms[i] = closestAtom
    # If one bead has only one heavy atom, include one more
    for i in atomPartitioning.values():
      if sum(x == i for x in atomPartitioning.values()) == 1:
        # Find bead
        lonelyBead = i
        # Voronoi to find closest atom
        closestBead = -1
        closestBeadDist = 10000.0
        for j in range(len(heavyAtomCoords)):
          if atomPartitioning[j] != lonelyBead:
            distBeadAt = np.linalg.norm(cgBeadCoords[lonelyBead] -
              heavyAtomCoords[j])
            # Only consider if it's closer, not a CG bead itself, and
            # the CG bead it belongs to has more than one other atom.
            if distBeadAt < closestBeadDist and \
              j != closestAtoms[atomPartitioning[j]] and \
              sum(x == atomPartitioning[j] for x in atomPartitioning.values()) > 2:
              closestBead = j
              closestBeadDist = distBeadAt
        if closestBead == -1:
          print("Error. Can't find an atom close to atom",lonelyBead)
          exit(1)
        atomPartitioning[closestBead] = lonelyBead
  return atomPartitioning

def checkAdditivity(beadTypes, mol, args):
  '''Check additivity assumption between sum of free energies of CG beads
  and free energy of whole molecule'''
  # If there's only one bead, don't check.
  if len(beadTypes) == 1:
    return True
  sumFrag = 0.0
  rings = False
  for bead in beadTypes:
    if bead[0] == "S": # bead belongs to a ring
      bead = bead[1:]
      rings = True
    sumFrag += topology.deltaFTypes[bead]
  # Wildman-Crippen logP
  wcLogP = rdMolDescriptors.CalcCrippenDescriptors(mol)[0]

  # Write out SMILES string of entire molecule
  tmpfile = 'tmp-auto-martini.smi'
  sw = Chem.rdmolfiles.SmilesWriter(tmpfile, nameHeader='')
  Chem.rdmolfiles.SmilesWriter.write(sw, mol)
  sw.close()
  # Read file
  try:
    f = open(tmpfile,'r')
    s = f.readlines()
    f.close()
  except IOError as e:
    print("Error. Can't read file", tmpfile)
    print(e)
    exit(1)
  os.remove(tmpfile)
  wholeMoldG = topology.smi2alogps(s[1].split()[0], wcLogP, "MOL", args, True)
  mad = math.fabs((wholeMoldG - sumFrag)/wholeMoldG)

  if args.verbose:
    print("; Mapping additivity assumption ratio: {0:7.4f} ({1:7.4f} vs {2:7.4f})".format(mad,wholeMoldG,sumFrag))
  if (not rings and mad < 0.5) or (rings):
    return True
  else:
    return False

if __name__ == '__main__':

# Parse command-line options
  parser = argparse.ArgumentParser(description='Map atomistic structure to MARTINI mapping', epilog='Tristan BEREAU (2014)')
  group = parser.add_mutually_exclusive_group(required=True)
  group.add_argument('--sdf', dest='sdf', type=str, required=False,
        help='SDF file of atomistic coordinates')
  group.add_argument('--smi', dest='smi', type=str, required=False,
        help='SMILES string of atomistic structure')
  parser.add_argument('--mol', dest='molname', type=str, required=True,
        help='Name of CG molecule')
  parser.add_argument('--xyz', dest='xyz', type=str,
        help='output atomistic structure to OUT.xyz file')
  parser.add_argument('--gro', dest='gro', type=str,
        help='output CG structure to OUT.gro file')
  parser.add_argument('--map', dest='map', type=str,
        help='output details of mapping')
  parser.add_argument('--verbose', dest='verbose', action='store_true',
        help='verbose')
  parser.add_argument('--fpred', dest='forcepred', action='store_true',
        help='verbose')

  args = parser.parse_args()

  if args.sdf:
    # Generate molecule's structure from SDF
    mol = topology.genMoleculeSDF(args.sdf)
  else:
    mol = topology.genMoleculeSMI(args.smi)

  feats = extractFeatures(mol)
  # Identify ring-type atoms
  ringAtoms = getRingAtoms(feats, args)
  # Get Hbond information
  hbondA = getHbondA(feats)
  hbondD = getHbondD(feats)

  # for feat in feats:
  #   print feat.GetFamily(),feat.GetType(),feat.GetAtomIds()

  # Flatten list of ring atoms
  ringAtomsFlat = list(chain.from_iterable(ringAtoms))

  # Optimize CG bead positions -- keep all possibilities in case something goes
  # wrong later in the code.
  listCGBeads,listBeadPos = findBeadPos(mol, feats, ringAtoms, args)

  # Loop through best 1% cgBeads and avgPos
  maxAttempts = int(math.ceil(0.5*len(listCGBeads)))

  if args.verbose:
    print("; Max. number of attempts:",maxAttempts)
  attempt = 0
  while attempt < maxAttempts:
    cgBeads = listCGBeads[attempt]
    beadPos  = listBeadPos[attempt]
    success  = True
    # Extract atom coordinates of heavy atoms
    heavyAtomCoords = getHeavyAtomCoords(mol) # remains constant
    # Extract position of CG beads
    cgBeadCoords = getCGBeadCoords(mol, cgBeads, beadPos, ringAtomsFlat) # changes with each iteration

    # Partition atoms into CG beads
    atomPartitioning = voronoiAtoms(cgBeadCoords, heavyAtomCoords)

    if args.verbose:
      print("; Atom partitioning:", atomPartitioning)

    atomNames, beadTypes = topology.printAtoms(cgBeads, mol, atomPartitioning, ringAtoms,
      ringAtomsFlat, hbondA, hbondD, args, True)

    if atomNames == []:
      success = False

    # Check additivity between fragments and entire molecule
    if not checkAdditivity(beadTypes, mol, args):
      success = False
    # Bond list
    #try:
    bondList, constList = topology.printBonds(cgBeads, mol, cgBeadCoords, ringAtoms, atomPartitioning, True) # this is so stupid ... how do we know if pringBonds technically fails?!
    #except:
    #  success = False

    if success:
      topology.printHeader(args)
      atomNames, beadTypes = topology.printAtoms(cgBeads, mol, atomPartitioning, ringAtoms,
        ringAtomsFlat, hbondA, hbondD, args, False)
      bondList, constList = topology.printBonds(cgBeads, mol, cgBeadCoords, ringAtoms, atomPartitioning, False)
      topology.printAngles(cgBeads, mol, atomPartitioning, cgBeadCoords,
        bondList, constList, ringAtoms)
      topology.printDihedrals(cgBeads,bondList,constList,ringAtoms,cgBeadCoords)
      # We've reached all the way here, exit the while loop
      attempt = maxAttempts+1
    else:
      attempt += 1
      #print("Running Iteration {}/{}".format(attempt, maxAttempts))

  if attempt == maxAttempts:
    err = "; ERROR: no successful mapping found.\n" + \
      "; Try running with the '--fpred' and/or '--verbose' options.\n"
    sys.stderr.write(err)
    sys.exit(0)

  # Optional atomistic output to XYZ file
  if args.xyz:
    output.outputXYZ(mol, args)
  # Optional CG output to GRO file
  if args.gro:
    output.outputGRO(cgBeadCoords, atomNames, args)
  # Optional output of mapping
  if args.map:
    output.outputMap(atomPartitioning, args)

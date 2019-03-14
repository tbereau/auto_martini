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
  - function printAtoms() uses undefined objects "hbondA" and "hbondD"
  - function checkAdditivity() defined "mad" variable which is already defined as a function
  - function genMoleculeSDF() did not use "sdf" input argument
  - function printBonds() used undefined "atomPartitioning" object

  TODO: make this run in Python 3

'''

import sys,os
import argparse
import math
import numpy as np
from rdkit import Chem
from rdkit.Chem import ChemicalFeatures
from rdkit.Chem import rdchem
from rdkit.Chem import rdMolDescriptors
from rdkit import RDConfig
from itertools import chain
import lxml

import cProfile

import output
import topology
from optimization import findBeadPos, getHeavyAtomCoords, voronoiAtoms

# For feature extraction
fdefName = os.path.join(RDConfig.RDDataDir,'BaseFeatures.fdef')
factory = ChemicalFeatures.BuildFeatureFactory(fdefName)


def getHbondA(feats):
  '''Get Hbond acceptor information'''
  hbondA = []
  for feat in feats:
    if feat.GetFamily() == "Acceptor":
      for i in feat.GetAtomIds():
        if i not in hbondA:
          hbondA.append(i)
  return hbondA

def getHbondD(feats):
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
    if feat.GetType() in ["RH6_6","RH5_5","RH4_4","RH3_3","Arom5","Arom6","Arom7","Arom8"]:
      newRing = []
      for at in feat.GetAtomIds():
          newRing.append(at)
      if newRing not in ringAtoms:
        ringAtoms.append(newRing)
  if args.verbose:
    print("; ring atoms:",ringAtoms)
  return ringAtoms

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

def run(mol, feats, ringAtoms, args):

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

    print("Running Iteration {}/{}".format(attempt, maxAttempts))

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
    try:
      bondList, constList = topology.printBonds(cgBeads, mol, cgBeadCoords, ringAtoms, atomPartitioning, True) # this is so stupid ... how do we know if pringBonds technically fails?!
    except:
      success = False

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

  run(mol, feats, ringAtoms, args)

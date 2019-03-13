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
'''

from rdkit import Chem
from rdkit.Chem import rdchem
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem import AllChem

import os, sys
from collections import defaultdict
from bs4 import BeautifulSoup
import requests
import math
import numpy as np

# Measured octanol/water free energies from MARTINI.
# Data used TI, not partitioning of Marrink et al. JPCB 2007.
deltaFTypes = {
  'Qda': -15.04,
  'Qa' : -15.04,
  'Qd' : -15.04,
  'Q0' : -22.35,
  'P5' :  -8.88,
  'P4' :  -9.30,
  'P3' :  -8.81,
  'P2' :  -3.85,
  'P1' :  -2.26,
  'Nda':   2.49,
  'Na' :   2.49,
  'Nd' :   2.49,
  'N0' :   4.22,
  'C5' :   6.93,
  'C4' :  10.14,
  'C3' :  12.26,
  'C2' :  13.74,
  'C1' :  14.20,
}

def genMoleculeSMI(smi):
  '''Generate mol object from smiles string'''
  smi = smi.upper()
  mol = Chem.MolFromSmiles(smi)
  mol = Chem.AddHs(mol)
  AllChem.EmbedMolecule(mol, randomSeed = 1, useRandomCoords=True) # Set Seed for random coordinate generation = 1.
  try:
    AllChem.UFFOptimizeMolecule(mol)
  except ValueError:
    print("Supply input args")
  return mol

def genMoleculeSDF(sdf):
  '''Generate mol object from SD file'''
  suppl = Chem.SDMolSupplier(sdf)
  if len(suppl) > 1:
    print("Error. Only one molecule may be provided.")
    sys.exit(1)
  for mol in suppl:
    if mol is None:
      print("Error. Can't read molecule.")
      sys.exit(1)
  return mol

def getCharge(mol):
  '''Get net charge of molecule'''
  return Chem.rdmolops.GetFormalCharge(mol)

def smi2alogps(smi, wcLogP, bead, args, trial=False):
  '''Returns water/octanol partitioning free energy
  according to ALOGPS'''

  session = requests.session()
  req = session.get('http://vcclab.org/web/alogps/calc?SMILES=' + str(smi))
  doc = BeautifulSoup(req.content, "lxml")
  soup = doc.prettify()
  foundMol1 = False

  for line in soup.split("\n"):
    line = line.split()
    if "mol_1" in line:
      logP = float(line[line.index('mol_1')+1])
      foundMol1 = True
      break
      
  if not foundMol1:
    # If we're forcing a prediction, use Wildman-Crippen
    if args.forcepred:
      if trial == True:
        wrn = "; Warning: bead ID " + str(bead) + \
          " predicted from Wildman-Crippen. Fragment " + str(smi) + "\n"
        sys.stderr.write(wrn)
      logP = wcLogP
    else:
      if args.verbose:
        print("ALOGPS can't predict fragment:", smi)
      exit(1)

  return convertLogK(logP)

def convertLogK(logK):
  '''Convert log_{10}K to free energy (in kJ/mol)'''
  return 0.008314*300.0*logK/math.log10(math.exp(1))

def mad(Type, deltaF):
  '''Mean absolute difference between type Type and deltaF'''
  return math.fabs(deltaFTypes[Type] - deltaF)

def determineBeadType(deltaF, charge, hbondA, hbondD, inRing):
  '''Determine CG bead type from deltaF value, charge,
  and hbond acceptor, and donor'''
  # if args.verbose:
  #   print "; dF:",deltaF/4.2,'kcal/mol'
  if charge < -1 or charge > +1:
    print("Charge is too large:", charge)
    print("No adequate force-field parameter")
    exit(1)
  beadType = []
  if inRing:
    # We're in a ring, divide free energy by 3
    # (average number of beads per ring)
    if abs(deltaF) > 0.1:
      deltaF *= 2/3.
  error = 0.0
  if charge != 0:
    # The compound has a +/- charge -> Q type
    error = mad('Qda',deltaF)
    if hbondA > 0 and hbondD > 0:
      beadType = "Qda"
    elif hbondA > 0 and hbondD == 0:
      beadType = "Qa"
    elif hbondA == 0 and hbondD > 0:
      beadType = "Qd"
    else:
      beadType = "Q0"
  else:
    # Neutral group
    # Use Hbond information only if we're close to Nda, Na, Nd types
    error = mad('Nda',deltaF)
    if error < 3.0 and (hbondA > 0 or hbondD > 0):
      if hbondA > 0 and hbondD > 0:
        beadType = "Nda"
      elif hbondA > 0 and hbondD == 0:
        beadType = "Na"
      elif hbondA == 0 and hbondD > 0:
        beadType = "Nd"
    else:
      # all other cases. Simply find the atom type that's closest in
      # free energy.
      otherTypes = ['P5','P4','P3','P2','P1','N0','C5','C4','C3','C2','C1']
      minError = 1000.0
      for cgType in otherTypes:
        tmpError = mad(cgType,deltaF)
        if tmpError < minError:
          minError = tmpError
          beadType = cgType
      error = minError
  # if error > 5:
  #   print "Warning: large error between beead type and logK value:"
  #   print " Type {:s} ({:5.2f}) vs. {:5.2f}".format(beadType,
  #     deltaFTypes[beadType],deltaF)
  if inRing:
    beadType = "S" + beadType
  return beadType

def printHeader(args):
  print(";;;; GENERATED WITH auto-martini")
  if args.smi:
    print("; INPUT SMILES:",args.smi)
  else:
    print(";",args.sdf)
  print("; Tristan Bereau (2014)")
  print("")
  print('[moleculetype]')
  print('; molname       nrexcl')
  print('  {:5s}         2'.format(args.molname))
  print('')
  print('[atoms]')
  print('; id    type    resnr   residu  atom    cgnr    charge  smiles')
  return

def printAtoms(cgBeads, mol, atomPartitioning, ringAtoms, ringAtomsFlat, hbondA, hbondD, args, trial=False):
  '''print CG Atoms in itp format'''
  atomNames = []
  beadTypes = []

  for bead in range(len(cgBeads)):
    # Determine SMI of substructure
    try:
      smiFrag, wcLogP, charge = substruct2smi(mol, atomPartitioning, bead, cgBeads, ringAtoms)
    except:
      raise

    atomName = ""

    for character, count in sorted(letterOccurrences(smiFrag).items()):
      try:
          float(character)
          pass
      except ValueError:
        if count == 1:
          atomName += "{:s}".format(character)
        else:
          atomName += "{:s}{:s}".format(character, str(count))
    # Get charge for smiFrag
    molFrag = genMoleculeSMI(smiFrag)
    chargeFrag = getCharge(molFrag)
    # Extract ALOGPS free energy
    try:
      if chargeFrag == 0:
        alogps = smi2alogps(smiFrag, wcLogP, bead+1, args, trial)
      else:
        alogps = 0.0
    except:
      raise
    hbondAFlag = 0
    for at in hbondA:
      if atomPartitioning[at] == bead:
        hbondAFlag = 1
        break
    hbondDFlag = 0
    for at in hbondD:
      if atomPartitioning[at] == bead:
        hbondDFlag = 1
        break
    inRing = True if cgBeads[bead] in ringAtomsFlat else False
    beadType = determineBeadType(alogps, charge, hbondAFlag, hbondDFlag, inRing)
    atomName = ""
    nameIndex = 0

    while atomName in atomNames or nameIndex == 0:
      nameIndex += 1
      atomName = "{:1s}{:02d}".format(beadType[0],nameIndex)
    atomNames.append(atomName)

    if trial == False:
      print('  {:<5d} {:5s}   1       {:5s}   {:7s} {:<5d}  {:2d}     ; {:s}'.format(bead+1, beadType, args.molname, atomName, bead+1, charge, smiFrag))

    beadTypes.append(beadType)

  return atomNames, beadTypes


def printBonds(cgBeads, mol, cgBeadCoords, ringAtoms, atomPartitioning, trial=False):
  '''print CG bonds in itp format'''
  if trial == False:
    print("")
  # Bond information
  bondList = []
  constList = []
  if len(cgBeads) > 1:
    for i in range(len(cgBeads)):
      for j in range(i+1,len(cgBeads)):
        dist = np.linalg.norm(cgBeadCoords[i]-cgBeadCoords[j])*0.1
        if dist < 0.65:
          # Are atoms part of the same ring
          inRing = False
          for ring in ringAtoms:
            if cgBeads[i] in ring and cgBeads[j] in ring:
              inRing = True
              break
          if inRing:
            constList.append([i,j,dist])
          else:
            # Check that the bond is not too short
            if dist < 0.2:
              raise NameError('Bond too short')
            # Look for a bond between an atom of i and an atom of j
            foundConnection = False
            atomsInBeadI = []
            for aa in atomPartitioning.keys():
              if atomPartitioning[aa] == i:
                atomsInBeadI.append(aa)
            atomsInBeadJ = []
            for aa in atomPartitioning.keys():
              if atomPartitioning[aa] == j:
                atomsInBeadJ.append(aa)
            for ib in range(len(mol.GetBonds())):
              abond = mol.GetBondWithIdx(ib)
              if (abond.GetBeginAtomIdx() in atomsInBeadI and \
                abond.GetEndAtomIdx() in atomsInBeadJ) or \
                (abond.GetBeginAtomIdx() in atomsInBeadJ and \
                abond.GetEndAtomIdx() in atomsInBeadI):
                foundConnection = True
            if foundConnection:
              bondList.append([i,j,dist])

    for ring in ringAtoms:
      # Only keep one bond between a ring and a given external bead
      for i in range(len(cgBeads)):
        at = cgBeads[i]
        if at not in ring:
          bondsToRing = []
          for b in bondList:
            if (cgBeads[b[0]] in ring and b[1] == at) or \
              (b[0] == at and cgBeads[b[1]] in ring):
              bondsToRing.append(b)
          # keep closest
          closestBond = [-1,-1,1000.0]
          for r in range(len(bondsToRing)):
            if bondsToRing[r][2] < closestBond[2]:
              closestBond = bondsToRing[r]
          # Delete the other bonds
          for b in bondsToRing:
            if b != closestBond:
              bondList.remove(b)
      beadsBondedToRing = []
      for i in range(len(cgBeads)):
        atomsInBead = []
        for key, val in atomPartitioning.iteritems():
          if val == i:
            atomsInBead.append(key)
        for b in bondList:
          if (b[0] in ring and b[1] in atomsInBead) or \
            (b[0] in atomsInBead and b[1] in ring):
            beadsBondedToRing.append(i)
      # Delete bond between 2 beads if they're both linked
      # to the same ring.
      for i in range(len(cgBeads)):
        for j in range(i+1,len(cgBeads)):
          if cgBeads[i] in beadsBondedToRing and \
            cgBeads[j] in beadsBondedToRing:
            for b in bondList:
              if (b[0] == i and b[1] == j) or \
                (b[0] == j and b[1] == i):
                bondList.remove(b)

    # Replace bond by constraint if both atoms have constraints
    # to the same third atom
    bondListIdx = 0
    while bondListIdx < len(bondList):
      b = bondList[bondListIdx]
      for i in range(len(cgBeads)):
        constI = False
        constJ = False
        for c in constList:
          if (c[0] == b[0] and c[1] == i) or \
            (c[0] == i and c[1] == b[0]):
            constI = True
          if (c[0] == b[1] and c[1] == i) or \
            (c[0] == i and c[1] == b[1]):
            constJ = True
        if constI and constJ:
          constList.append(b)
          bondList.remove(b)
      bondListIdx += 1

    # Go through list of constraints. If we find an extra
    # possible constraint between beads that have constraints,
    # add it.
    beadsWithConst = []
    for c in constList:
      if c[0] not in beadsWithConst:
        beadsWithConst.append(c[0])
      if c[1] not in beadsWithConst:
        beadsWithConst.append(c[1])
    beadsWithConst = sorted(beadsWithConst)
    for i in range(len(beadsWithConst)):
      for j in range(1+i,len(beadsWithConst)):
        constExists = False
        for c in constList:
          if (c[0] == i  and c[1] == j) or \
            (c[0] == j and c[1] == i):
            constExists = True
            break
        if not constExists:
          dist = np.linalg.norm(cgBeadCoords[i]-cgBeadCoords[j])*0.1
          if dist < 0.35:
            # Are atoms part of the same ring
            inRing = False
            for ring in ringAtoms:
              if cgBeads[i] in ring and cgBeads[j] in ring:
                inRing = True
                break
            # Check that it's not in the bond list
            inBondList = False
            for b in bondList:
              if (b[0] == i and b[1] == j) or \
                (b[0] == j and b[0] == i):
                inBondList = True
                break
            if not inBondList and inRing:
              constList.append([i,j,dist])

    if trial == False:
      if len(bondList) > 0:
        print("[bonds]")
        print("; i j   funct   length  force.c.")
        for b in bondList:
          # Make sure atoms in bond are not part of the same ring
          print("  {:d} {:d}   1       {:4.2f}    1250".format(b[0]+1,b[1]+1,b[2]))
        print("")
      if len(constList) > 0:
        print("[constraints]")
        print(";  i   j     funct   length")
        for c in constList:
          print("   {:<3d} {:<3d}   1       {:4.2f}".format(c[0]+1,c[1]+1,c[2]))
        print("")
      # Make sure there's at least a bond to every atom
      for i in range(len(cgBeads)):
        bondToI = False
        for b in bondList + constList:
          if i in [b[0],b[1]]:
            bondToI = True
        if not bondToI:
          print("Error. No bond to atom",i+1)
          exit(1)
  return bondList, constList

def printAngles(cgBeads, mol, atomPartitioning, cgBeadCoords,
  bondList, constList, ringAtoms):
  '''print CG angles in itp format'''
  if len(cgBeads) > 2:
    # Angles
    angleList = []
    for i in range(len(cgBeads)):
      for j in range(len(cgBeads)):
        for k in range(len(cgBeads)):
          # Only up to 2 atoms can be ring-like
          allInRing = False
          for ring in ringAtoms:
            if cgBeads[i] in ring and cgBeads[j] in ring and \
              cgBeads[k] in ring:
              allInRing = True
              break
          # Forbid all atoms linked by constraints
          allConstraints = False
          ijBonded = False
          jkBonded = False
          ijConst  = False
          jkConst  = False
          for b in bondList + constList:
            if i in [b[0],b[1]] and j in [b[0],b[1]]:
              ijBonded = True
              if b in constList:
                ijConst = True
            if j in [b[0],b[1]] and k in [b[0],b[1]]:
              jkBonded = True
              if b in constList:
                jkConst = True
          if ijConst and jkConst:
            allConstraints = True
          if not allInRing and ijBonded and jkBonded and \
            i != j and j != k and i != k and \
            not allConstraints:
            # Measure angle between i, j, and k.
            angle = 180./math.pi*math.acos(
              np.dot(cgBeadCoords[i]-cgBeadCoords[j],
              cgBeadCoords[k]-cgBeadCoords[j])/(
              np.linalg.norm(cgBeadCoords[i]-cgBeadCoords[j]) *
              np.linalg.norm(cgBeadCoords[k]-cgBeadCoords[j])))
            # Look for any double bond between atoms belonging to these CG beads.
            atomsInFragment = []
            for aa in atomPartitioning.keys():
              if atomPartitioning[aa] == j:
                atomsInFragment.append(aa)
            forcConst = 25.0
            for ib in range(len(mol.GetBonds())):
              abond = mol.GetBondWithIdx(ib)
              if abond.GetBeginAtomIdx() in atomsInFragment and \
                abond.GetEndAtomIdx() in atomsInFragment:
                bondType = mol.GetBondBetweenAtoms(abond.GetBeginAtomIdx(),
                  abond.GetEndAtomIdx()).GetBondType()
                if bondType == Chem.rdchem.BondType.DOUBLE:
                  forcConst = 45.0
            newAngle = True
            for a in angleList:
              if i in a and j in a and k in a:
                newAngle = False
            if newAngle:
              angleList.append([i,j,k,angle,forcConst])
    if len(angleList) > 0:
      print("[angles]")
      print("; i j k         funct   angle   force.c.")
      for a in angleList:
        print("  {:d} {:d} {:d}         2       {:<5.1f}  {:5.1f}".format(a[0]+1,a[1]+1,a[2]+1,a[3],a[4]))
      print("")
  return

def printDihedrals(cgBeads,bondList,constList,ringAtoms,cgBeadCoords):
  '''Print CG dihedrals in itp format'''
  if len(cgBeads) > 3:
    # Dihedrals
    dihedList = []
    # Three ring atoms and one non ring
    for i in range(len(cgBeads)):
      for j in range(len(cgBeads)):
        for k in range(len(cgBeads)):
          for l in range(len(cgBeads)):
            if i != j and i != k and i != l \
              and j != k and j != l and k != l:
              # 3 atoms need to be ring like (in one ring!)
              threeInRing = False
              for ring in ringAtoms:
                if [[cgBeads[i] in ring], [cgBeads[j] in ring],
                  [cgBeads[k] in ring], [cgBeads[l] in ring]].count([True]) >= 3:
                  threeInRing = True
                  break
              for b in constList:
                if i in [b[0],b[1]] and j in [b[0],b[1]]:
                  ijConst = True
                if j in [b[0],b[1]] and k in [b[0],b[1]]:
                  jkConst = True
                if k in [b[0],b[1]] and l in [b[0],b[1]]:
                  klConst = True
              # Distance criterion--beads can't be far apart
              disThres = 0.35
              closeEnough = False
              if np.linalg.norm(cgBeadCoords[i]-cgBeadCoords[j])*0.1 < disThres \
                and np.linalg.norm(cgBeadCoords[j]-cgBeadCoords[k])*0.1 < disThres \
                and np.linalg.norm(cgBeadCoords[k]-cgBeadCoords[l])*0.1 < disThres:
                closeEnough = True
              alreadyDih = False
              for dih in dihedList:
                if dih[0] == l and dih[1] == k and \
                  dih[2] == j and dih[3] == i:
                  alreadyDih = True
                  break
              if threeInRing and closeEnough and not alreadyDih:
                r1 = cgBeadCoords[j]-cgBeadCoords[i]
                r2 = cgBeadCoords[k]-cgBeadCoords[j]
                r3 = cgBeadCoords[l]-cgBeadCoords[k]
                p1 = np.cross(r1, r2)/(np.linalg.norm(r1)*
                  np.linalg.norm(r2))
                p2 = np.cross(r2, r3)/(np.linalg.norm(r2)*
                  np.linalg.norm(r3))
                r2 = r2 / np.linalg.norm(r2)
                cosphi = np.dot(p1,p2)
                sinphi = np.dot(r2,np.cross(p1,p2))
                angle = 180./math.pi*np.arctan2(sinphi,cosphi)
                forcConst = 10.0
                dihedList.append([i,j,k,l,angle,forcConst])
    if len(dihedList) > 0:
      print("[dihedrals]")
      print(";  i     j    k    l   funct   angle  force.c.")
      for d in dihedList:
        print("   {:d}     {:d}    {:d}    {:d}       2     {:<5.1f}  {:5.1f}".format(d[0]+1,d[1]+1,d[2]+1,d[3]+1,d[4],d[5]))
      print("")
  return

def substruct2smi(mol, atomPartitioning, cgBead, cgBeads, ringAtoms):
  '''Substructure to smiles conversion; also output Wildman-Crippen logP;
     and charge of group.'''
  frag = rdchem.EditableMol(mol)
  numAtoms = mol.GetConformer().GetNumAtoms()
  # First delete all hydrogens
  for i in range(numAtoms):
    if mol.GetAtomWithIdx(i).GetSymbol() == "H":
      # find atom from coordinates
      submol = frag.GetMol()
      for j in range(submol.GetConformer().GetNumAtoms()):
        if mol.GetConformer().GetAtomPosition(i)[0] == \
          submol.GetConformer().GetAtomPosition(j)[0]:
          frag.RemoveAtom(j)
  # Identify atoms involved in same ring as cgBead (only one ring)
  atomsInRing = []
  for ring in ringAtoms:
    if cgBeads[cgBead] in ring:
      atomsInRing = ring[:] ;# CHANGED
      break
  # Add atoms off the ring that belong to the fragment.
  for atom in atomsInRing:
    if atom == cgBeads[cgBead]:
      for atp in atomPartitioning.keys():
        if atomPartitioning[atp] == atomPartitioning[atom] and atp not in atomsInRing:
          atomsInRing.append(atp)
  # Then heavy atoms that aren't part of the CG bead (except those
  # involved in the same ring).
  for i in atomPartitioning.keys():
    submol = frag.GetMol()
    if atomPartitioning[i] != cgBead and i not in atomsInRing:
      # find atom from coordinates
      submol = frag.GetMol()
      for j in range(submol.GetConformer().GetNumAtoms()):
        if mol.GetConformer().GetAtomPosition(i)[0] == \
          submol.GetConformer().GetAtomPosition(j)[0]:
          frag.RemoveAtom(j)
  # Wildman-Crippen logP
  wcLogP = rdMolDescriptors.CalcCrippenDescriptors(frag.GetMol())[0]
  tmpfile = 'tmp-auto-martini.smi'
  sw = Chem.rdmolfiles.SmilesWriter(tmpfile,nameHeader='')
  # Charge -- look at atoms that are only part of the bead (no ring rule)
  chg = 0
  for i in atomPartitioning.keys():
    if atomPartitioning[i] == cgBead:
      chg += mol.GetAtomWithIdx(i).GetFormalCharge()
  # Chem.rdmolfiles.SmilesWriter.write(sw,frag.GetMol())
  Chem.rdmolfiles.SmilesWriter.write(sw,
    Chem.rdmolops.AddHs(frag.GetMol(), addCoords=True )) # Andrew: problem is this seems to write the SMILES code in lower-case letters
  sw.close()
  # Read file
  try:
    f = open(tmpfile,'r')
    s = f.readlines()
    f.close()
  except IOError as e:
    print("Error. Can't read file",tmpfile)
    print(e)
    exit(1)
  os.remove(tmpfile)

  smiFrag, wcLogP, charge = s[1].split()[0], wcLogP, chg

  return smiFrag, wcLogP, charge 

def letterOccurrences(string):
  frequencies = defaultdict(lambda: 0)
  for character in string:
    if character.isalnum():
      frequencies[character.upper()] += 1
  return frequencies
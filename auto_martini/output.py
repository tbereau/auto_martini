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

import os

def outputXYZ(mol, args):
  '''Ouput XYZ file of molecule object'''
  numAtoms = mol.GetConformer().GetNumAtoms()
  if args.xyz[-4:] != ".xyz":
    args.xyz = args.xyz + ".xyz"
  try:
    with open(args.xyz,'w') as f:
      f.write(str(numAtoms) + "\n")
      f.write(" " + args.xyz[:-4] + "\n")
      for i in range(numAtoms):
        f.write("{:2s}  {:7.4f} {:7.4f} {:7.4f}\n".format(
          mol.GetAtomWithIdx(i).GetSymbol(),
          mol.GetConformer().GetAtomPosition(i)[0],
          mol.GetConformer().GetAtomPosition(i)[1],
          mol.GetConformer().GetAtomPosition(i)[2]))
      f.write("\n")
      f.close()
  except IOError:
    print("Can't write to file " + args.xyz)
    exit(1)
  return

def outputGRO(beads, beadNames, args):
  '''Output GRO file of CG structure'''
  numBeads = len(beads)
  if len(beads) != len(beadNames):
    print("Error. Incompatible number of beads and bead names.")
    exit(1)
  if args.gro[-4:] != ".gro":
    args.gro = args.gro + ".gro"
  try:
    with open(args.gro,'w') as f:
      f.write("{:s} generated from {:s}\n".format(
        args.molname,os.path.basename(__file__)))
      f.write("{:5d}\n".format(numBeads))
      for i in range(numBeads):
        f.write("{:5d}{:<6s} {:3s}{:5d}{:8.3f}{:8.3f}{:8.3f}\n".format(
          i+1, args.molname, beadNames[i], i+1, beads[i][0]/10.,
          beads[i][1]/10., beads[i][2]/10.))
      f.write("{:10.5f}{:10.5f}{:10.5f}\n".format(10.,10.,10.))
      f.close()
  except IOError:
    print("Can't write to file " + args.gro)
    exit(1)
  return

def outputPDB(mol, cgBeads, args):
  '''Output PDB file of AA/CG structure'''
  for i in range(len(cgBeads)):
    ati = mol.GetAtomWithIdx(i)
    pdbres = Chem.rdchem.AtomPDBResidueInfo(ati)
    print(pdbres.GetResidueName())
    print(ati.GetSmarts())
  pw = Chem.rdmolfiles.PDBWriter(args.pdb)
  Chem.rdmolfiles.PDBWriter.write(pw, mol)
  pw.close()
  return

def outputMap(atomPartitioning, args):
  '''Output the all-atom mapping for each bead'''
  beadMap = {}
  # Turn around the atomPartitioning so it is CG -> AA, rather than AA -> CG
  for atom, bead in atomPartitioning.iteritems() :
    if bead not in beadMap:
      beadMap[bead] = []
    beadMap[bead].append(atom)
  try:
    with open(args.map, 'w') as f:
      for bead in sorted(beadMap.keys()):
        f.write("{:d} : ".format(bead+1))
        f.write(" ".join(["{:d}".format(a+1) for a in beadMap[bead]]))
        f.write("\n")
      f.close()
  except IOError:
    print("Can't write to file " + args.map)
    exit(1)
  return
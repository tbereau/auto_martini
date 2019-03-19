# !/usr/bin/python
#  -*- coding: utf8 -*-

'''
  Created on March 17, 2019 by Andrew Abi-Mansour

          _    _ _______ ____    __  __          _____ _______ _____ _   _ _____ 
     /\  | |  | |__   __/ __ \  |  \/  |   /\   |  __ \__   __|_   _| \ | |_   _|
    /  \ | |  | |  | | | |  | | | \  / |  /  \  | |__) | | |    | | |  \| | | |  
   / /\ \| |  | |  | | | |  | | | |\/| | / /\ \ |  _  /  | |    | | | . ` | | |  
  / ____ \ |__| |  | | | |__| | | |  | |/ ____ \| | \ \  | |   _| |_| |\  |_| |_ 
 /_/    \_\____/   |_|  \____/  |_|  |_/_/    \_\_|  \_\ |_|  |_____|_| \_|_____|
                                                                                 
                                                                                 
  @originally written by: Tristan BEREAU (bereau at mpip-mainz.mpg.de)
  @modified by: Andrew Abi-Mansour (andrew.gaam at gmail.com)

  Auto_Martini: a tool for automatic MARTINI mapping and parametrization of small organic molecules                                               
    
  Auto_Martini is open-source, distributed under the terms of the GNU Public
  License, version 2 or later. It is distributed in the hope that it will
  be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
  of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. You should have
  received a copy of the GNU General Public License along with PyGran.
  If not, see http://www.gnu.org/licenses . See also top-level README
  and LICENSE files.

  Github link to original repo: https://github.com/tbereau/auto_martini
  '''

import argparse
import logging

from auto_martini.engine.solver import cg_molecule
from auto_martini.engine.topology import gen_molecule_smi, gen_molecule_sdf

import sys


def checkArgs(args):

  if not args.sdf and not args.smi:
    parser.error("run requires --sdf or --smi")
  
  if not args.molname:
    parser.error("run requires --mol")
  
  if not args.topfname:
    parser.error("run requires --top")

parser = argparse.ArgumentParser(prog='auto_martini', description='Generates Martini force field for atomistic structures of small organic molecules',
                                formatter_class=argparse.RawDescriptionHelpFormatter,
                                 epilog='''Developers:\n===========\nTristan Bereau (bereau [at] mpip-mainz.mpg.de)\nAndrew Abi-Mansour (andrew.gaam [at] gmail.com)''')

parser.add_argument('mode', type=str, help='run or test auto_martini', choices=['run','test'])

group = parser.add_mutually_exclusive_group(required=False)

group.add_argument('--sdf', dest='sdf', type=str, required=False, help='SDF file of atomistic coordinates')
group.add_argument('--smi', dest='smi', type=str, required=False, help='SMILES string of atomistic structure')
parser.add_argument('--mol', dest='molname', type=str, required=False, help='Name of CG molecule')
parser.add_argument('--aa', dest='aa', type=str, help='filename of all-atom structure .gro file')
parser.add_argument('--cg', dest='cg', type=str, help='filename of coarse-grained structure .gro file')
parser.add_argument('--top', dest='topfname', type=str, help='filename of output topology file')
parser.add_argument('-v', '--verbose', dest='verbose', action='count', default=0, help='increase verbosity')
parser.add_argument('--fpred', dest='forcepred', action='store_true', help='verbose')

if len(sys.argv)==1:
    parser.print_help(sys.stderr)
    sys.exit(1)

args = parser.parse_args()

if args.mode == 'run':
  checkArgs(args)

  if args.verbose >= 2:
      level = logging.DEBUG
  elif args.verbose >= 1:
      level = logging.INFO
  else:
      level = logging.WARNING

  # Generate molecule's structure from SDF or SMILES
  if args.sdf:
      mol = gen_molecule_sdf(args.sdf)
  else:
      mol = gen_molecule_smi(args.smi) 

  cg_molecule(mol, args.molname, args.topfname, args.aa, args.cg, args.forcepred)

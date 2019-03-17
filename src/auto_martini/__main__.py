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

  Auto_MARTINI: a tool for automatic MARTINI mapping and parametrization of small organic molecules                                               
    
  Auto_MARTINI is open-source, distributed under the terms of the GNU Public
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

from auto_martini.solver import cg_molecule
from auto_martini.topology import gen_molecule_sdf

parser = argparse.ArgumentParser(description='Map atomistic structure to MARTINI mapping',
                                 epilog='Tristan BEREAU (2014); Andrew Abi-Mansour (2019)')
group = parser.add_mutually_exclusive_group(required=True)
group.add_argument('--sdf', dest='sdf', type=str, required=False, help='SDF file of atomistic coordinates')
group.add_argument('--smi', dest='smi', type=str, required=False, help='SMILES string of atomistic structure')
parser.add_argument('--mol', dest='molname', type=str, required=True, help='Name of CG molecule')
parser.add_argument('--aa', dest='aa', type=str, help='output all-atom structure to .gro file')
parser.add_argument('--cg', dest='cg', type=str, help='output coarse-grained structure to .gro file')
parser.add_argument('-v', '--verbose', dest='verbose', action='count', default=0, help='increase verbosity')
parser.add_argument('--fpred', dest='forcepred', action='store_true', help='verbose')

args = parser.parse_args()
if args.verbose >= 2:
    level = logging.DEBUG
elif args.verbose >= 1:
    level = logging.INFO
else:
    level = logging.WARNING

logging.basicConfig(filename='auto_martini.log', format='%(asctime)s:%(levelname)s: %(message)s', level=level)

# Generate molecule's structure from SDF or SMILES
if args.sdf:
    mol = gen_molecule_sdf(args.sdf)
else:
    mol = gen_molecule_smi(args.smi) 

cg_molecule(mol, args.molname, args.aa, args.cg, args.forcepred)

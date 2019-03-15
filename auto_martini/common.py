# !/usr/bin/python
#  -*- coding: utf8 -*-

'''
  Created on March 14, 2019 by Andrew Abi-Mansour

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

from __future__ import print_function

import sys
import itertools
import os
import argparse
import math
import logging
import numpy as np
import six
import requests
from bs4 import BeautifulSoup
from collections import Counter
from rdkit import Chem
from rdkit import RDConfig
from rdkit.Chem import AllChem
from rdkit.Chem import ChemicalFeatures
from rdkit.Chem import rdchem
from rdkit.Chem import rdmolops
from rdkit.Chem import rdmolfiles
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem import rdMolTransforms
from itertools import chain
from collections import defaultdict
from operator import itemgetter
from sanifix4 import AdjustAromaticNs

import cProfile
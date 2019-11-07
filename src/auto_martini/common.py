'''
Created on March 14, 2019 by Andrew Abi-Mansour

This is the::

	     _   _   _ _____ ___    __  __    _    ____ _____ ___ _   _ ___ 
	    / \ | | | |_   _/ _ \  |  \/  |  / \  |  _ \_   _|_ _| \ | |_ _|
	   / _ \| | | | | || | | | | |\/| | / _ \ | |_) || |  | ||  \| || | 
	  / ___ \ |_| | | || |_| | | |  | |/ ___ \|  _ < | |  | || |\  || | 
	 /_/   \_\___/  |_| \___/  |_|  |_/_/   \_\_| \_\|_| |___|_| \_|___|                                                            
                                                                 
Tool for automatic MARTINI mapping and parametrization of small organic molecules

Developers::

	Tristan BEREAU (bereau at mpip-mainz.mpg.de)
	Kiran Kanekal (kanekal at mpip-mainz.mpg.de)
	Andrew Abi-Mansour (andrew.gaam at gmail.com)

AUTO_MARTINI is open-source, distributed under the terms of the GNU Public
License, version 2 or later. It is distributed in the hope that it will
be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. You should have
received a copy of the GNU General Public License along with PyGran.
If not, see http://www.gnu.org/licenses . See also top-level README
and LICENSE files.
'''

from __future__ import print_function
import logging

import sys
import itertools
import os
import math

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
from .sanifix4 import AdjustAromaticNs
from . import __version__

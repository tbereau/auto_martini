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

  BUGS found and fixed:

  - function substruct2smi() returned SMILES string in lower case letters
  - function printAtoms() uses undefined objects "hbondA" and "hbondD"
  - function checkAdditivity() defined "mad" variable which is already defined as a function
  - function genMoleculeSDF() did not use "sdf" input argument
  - function printBonds() used undefined "atomPartitioning" object
  - function Substructure() no longer creates then reads "tmp-auto-martini.smi" file


  TODO: make this run in Python 3

'''

from setuptools import setup, find_packages

from src.auto_martini import __version__

try:
	from Cython.Build import cythonize
	import numpy
	optimal_list = cythonize("src/auto_martini/optimization.pyx")
	include_dirs = [numpy.get_include()]
except:
	optimal_list = []
	include_dirs = []

setup(
    name = "auto_martini",
    version = __version__,
    author = ["Tristan Bereau", "Andrew Abi-Mansour"],
    author_email = ["bereau [at] mpip-mainz.mpg.de", "andrew.gaam [at] gmail.com"],
    description = ("A tool for automatic MARTINI mapping and parametrization of small organic molecules "),
    license = "GPL v2",
    keywords = "Coarse-grained Molecular Dynamics, MARTINI force field",
    url = "https://github.com/Andrew-AbiMansour/Auto_MARTINI",
    packages=['auto_martini'],
    package_dir={'auto_martini': 'src/auto_martini'},
    package_data={'auto_martini': ['examples/*.sdf', 'LICENSE', 'README.md']},
    include_package_data=True,
    install_requires=['numpy', 'bs4', 'pytool'],
    classifiers=[
        "Development Status :: 2 - Pre-Alpha",
        "Topic :: Utilities",
        "License :: OSI Approved :: GNU General Public License v2 (GPLv2)",
	"Programming Language :: Python :: 2.7",
	"Programming Language :: Python :: 3.6"
    ],
    zip_safe=False,
    ext_modules=optimal_list,
    include_dirs=include_dirs
)

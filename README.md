[![Build Status](https://circleci.com/gh/Andrew-AbiMansour/Auto_Martini.svg?branch=master)](https://circleci.com/gh/Andrew-AbiMansour/Auto_Martini)

Auto_MARTINI
============
*Developers:* 
Tristan Bereau (Max Planck Institute for Polymer Research, Mainz, Germany)  
Kiran Kanekal (Max Planck Institute for Polymer Research, Mainz, Germany) 
Andrew Abi-Mansour (Molecular Science Software Institute, Virginia Tech, Blacksburg, US) 

Automated MARTINI mapping and parametrization of small organic molecules. This fork restructures the original code to a modular and extensible form with numerous bug fixes and support for Python 3. For full documentation, click [here](https://andrew-abimansour.github.io/Auto_Martini/docs/html/index.html).

## Publication
For a detailed account of the software, see:

Bereau and Kremer, *J Chem Theory Comput*, DOI:10.1021/acs.jctc.5b00056 (2015) [link](http://dx.doi.org/10.1021/acs.jctc.5b00056)

[![DOI for Citing auto_martini](https://img.shields.io/badge/DOI-10.1021%2Facs.jctc.5b00056-blue.svg)](http://dx.doi.org/10.1021/acs.jctc.5b00056)

Please consider citing the paper if you find `auto_martini` useful in your research.

```
@article{bereau2015automartini,
author = {Bereau, Tristan and Kremer, Kurt},
title = {Automated parametrization of the coarse-grained MARTINI force field 
    for small organic molecules},
journal = {J Chem Theory Comput},
year = {2015},
volume = {11},
number = {6},
pages = {2783--2791},
doi = {10.1021/acs.jctc.5b00056}
}
```

## Installation & Testing
`auto-martini` is a python module that requires a number of dependencies:
* `numpy`: see http://docs.scipy.org/doc/numpy/user/install.html
* `rdkit`: see http://www.rdkit.org/docs/Install.html
* `bs4`: see http://www.crummy.com/software/BeautifulSoup/
* `requests`: see http://docs.python-requests.org/en/latest/user/install
* `lxml`: see https://github.com/lxml/lxml

For a python 2 installation, pip can be used to install all 5 depdendencies:
```
pip install numpy rdkit bs4 lxml requests
``` 
For optimal performance, however, we recommend you use python 3, which requires installing rdkit from its source code. Once all the dependencies are correctly installed, auto_martini can be run as a module via:
```
python -m auto_martini [mode] [options]
```
By default, mode is set to 'run', which computes the MARTINI forcefield for a given molecule. In order to make sure auto_martini runs correctly on your system, run:
```
python -m auto_martini test
```
To display the usage-information (help), either supply -h, --help, or nothing to auto_martini:
 
```
usage: auto_martini [-h] [--mode {run,test}] [--sdf SDF | --smi SMI]
                    [--mol MOLNAME] [--aa AA] [--cg CG] [--top TOPFNAME] [-v]
                    [--fpred]

Generates Martini force field for atomistic structures of small organic molecules

optional arguments:
  -h, --help         show this help message and exit
  --mode {run,test}  mode: run (compute FF) or test (validate)
  --sdf SDF          SDF file of atomistic coordinates
  --smi SMI          SMILES string of atomistic structure
  --mol MOLNAME      Name of CG molecule
  --aa AA            filename of all-atom structure .gro file
  --cg CG            filename of coarse-grained structure .gro file
  --top TOPFNAME     filename of output topology file
  -v, --verbose      increase verbosity
  --fpred            verbose

Developers:
===========
Tristan Bereau (bereau [at] mpip-mainz.mpg.de)
Kiran Kanekal (kanekal [at] mpip-mainz.mpg.de)
Andrew Abi-Mansour (andrew.gaam [at] gmail.com)
```

## Example
To coarse-grain a molecule, simply provide its SMILES code (option `--smi SMI`) or a .SDF file (option `'--sdf file.sdf`). You also need to provide a name for the CG molecule (not longer than 5 characters) using the `--mol` option.  For instance, to coarse grain [guanazole](http://pubchem.ncbi.nlm.nih.gov/summary/summary.cgi?cid=15078), you can either obtain/generate (e.g., from Open Babel) an SDF file:
```
python -m auto-martini --sdf guanazole.sdf --mol GUA --top GUA.itp
```
(the name GUA is arbitrary) or use its SMILES code within double quotes
```
python -m auto-martini --smi "N1=C(N)NN=C1N" --mol GUA --top GUA.itp
```
In case no problem arises, it will output the gromacs GUA.itp file:
```
;;;; GENERATED WITH auto-martini
; INPUT SMILES: N1=C(N)NN=C1N
; Tristan Bereau (2014)

[moleculetype]
; molname       nrexcl
  GUA           2

[atoms]
; id    type    resnr   residu  atom    cgnr    charge  smiles
  1     SP2     1       GUA     S01     1       0     ; Nc1ncnn1
  2     SP2     1       GUA     S02     2       0     ; Nc1ncnn1

[constraints]
;  i   j     funct   length
   1   2     1       0.21
```
Optionally, the code can also output a corresponding `.gro` file for the coarse-grained coordinates
```
python -m auto-martini --smi "N1=C(N)NN=C1N" --mol GUA --gro gua.gro --top GUA.itp
```
Atomistic coordinates can be output in XYZ format using the `--xyz output.xyz` option.

## Caveats

### Prediction algorithm

Since ALOGPS, the prediction algorithm for octanol/water partitioning, relies on whole fragments rather than individual atoms, the prediction of certain fragments can pose problem, e.g., small inorganic groups. In this case, `auto-martini` tries to parametrize alternative mappings. If none of them shows successful, the code will return an error.
```
; ERROR: no successful mapping found.
; Try running with the '--fpred' and/or '--verbose' options.
```
As mentioned in the error message, an alternative solution consists of relying on an atom-based partitioning coefficient prediction algorithm (Wildman-Crippen), which is less accurate but can predict any fragment.  In case the `--fpred` option is selected, only fragments for which ALOGPS fail will be predicted using Wildman-Crippen.

### Boost error

Some versions of Boost will fail to correctly exit at the end of the program, generating such output messages:
```
python: /usr/include/boost/thread/pthread/mutex.hpp:108: boost::mutex::~mutex(): Assertion `!posix::pthread_mutex_destroy(&m)' failed.
[1]    31433 abort (core dumped)  ./auto_martini --smi "N1=C(N)NN=C1N" --mol GUA
```
the results provided by the code are unaffected by this error message. Simply ignore it.

### RDKit outdated

Older RDKit versions will report the following error:
```
[...]
distBdAt = Chem.rdMolTransforms.GetBondLength(conf,i,beadId)
AttributeError: 'module' object has no attribute 'GetBondLength'
```
Simply update your version of RDKit. Most package managers will allow you to do this, unless you've installed RDKit from source.

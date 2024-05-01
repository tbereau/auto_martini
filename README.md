Auto_Martini
============

[![CI](https://github.com/tbereau/auto_martini/actions/workflows/CI.yaml/badge.svg)](https://github.com/tbereau/auto_martini/actions/workflows/CI.yaml)
[![codecov](https://codecov.io/gh/tbereau/auto_martini/branch/main/graph/badge.svg)](https://codecov.io/gh/tbereau/auto_martini/branch/main)

For reporting bugs or suggesting new features/improvements to the code, please open an [issue](https://github.com/tbereau/auto_martini/issues).

---

## What is Auto_Martini?
A toolkit that enables automatic generation of Martini forcefields for small organic molecules. 

For a detailed account of the software, see:

Bereau and Kremer, *J Chem Theory Comput*, DOI:10.1021/acs.jctc.5b00056 (2015)

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
pages = {2783-2791},
doi = {10.1021/acs.jctc.5b00056}
}
```

For full documentation, click [here](https://tbereau.github.io/auto_martini/docs/html/index.html).

## Developers 
* Tristan Bereau (University of Amsterdam, Netherlands)   
* Kiran Kanekal (Max Planck Institute for Polymer Research, Mainz, Germany)     
* Andrew Abi-Mansour

## Update to Python3
The `main` branch is now fully compatible with Python 3. For the original Python2-based version of the code used in the JCTC 2015 paper, see [branch original_jctc2015](https://github.com/tbereau/auto_martini/tree/original_jctc2015).

## Installation

### Installation with pip
You can install auto_martini with pip by running the following command from the source directory:
```bash
pip install .
```
If you do not wish to clone the repo, you can instead run:
```bash
pip install git+https://github.com/tbereau/auto_martini
```

### Installation with poetry
You can install auto_martini with [poetry](https://python-poetry.org) by running the following command from the source directory:
```bash
poetry install
```
This will install auto_martini in a virtual environment, which you can activate via:
```bash
poetry shell
```
To exit the environment, simply run `deactivate`.

### Installation with conda-lock
You can use [conda-lock](https://conda.github.io/conda-lock) to install auto_martini as well. From the source directory run:
```
conda-lock -f pyproject.toml  -k explicit --filename-template auto_martini-py3.11.conda.lock
```
This generates a conda lock file which you can use to create a new conda virtual environment:
```
conda create --name YOURENV --file auto_martini-py3.11.conda.lock 
```
`YOURENV` is the name of your conda env, and the lock filename can change depending on the python interpreter version.
Finally activate your conda envrionment via:
```
conda activate YOURENV
```
Now you can install the pkg with pip by running from the source directory:
```
pip install .
```
To exit the environment, simply run `conda deactivate`.

## Testing
To run the test cases and validate your installation, you will need to have [pytest](https://docs.pytest.org/en/stable/getting-started.html) 
installed. If you installed `auto_martini` with conda, then pytest should already be available in your environment.

To initiate testing, activate the virtual environment and run the following from the source directory:
```bash
pytest -v auto_martini/tests
```

All tests should pass within few minutes. If any of the tests fail, please open an [issue](https://github.com/tbereau/auto_martini/issues).

## Command-line Interface
You can invoke `auto_martini` from the command-line via:
```
python -m auto_martini [mode] [options]
```
By default, mode is set to 'run', which computes the MARTINI forcefield for a given molecule.

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
  --fpred            Atomic partitioning prediction

Developers:
===========
Tristan Bereau (bereau [at] mpip-mainz.mpg.de)
Kiran Kanekal (kanekal [at] mpip-mainz.mpg.de)
Andrew Abi-Mansour (andrew.gaam [at] gmail.com)
```

## Example
To coarse-grain a molecule, simply provide its SMILES code (option `--smi SMI`) or a .SDF file (option `'--sdf file.sdf`). You also need to provide a name for the CG molecule (not longer than 5 characters) using the `--mol` option.  For instance, to coarse grain [guanazole](http://pubchem.ncbi.nlm.nih.gov/summary/summary.cgi?cid=15078), you can either obtain/generate (e.g., from Open Babel) an SDF file:
```
python -m auto_martini --sdf guanazole.sdf --mol GUA --top GUA.itp
```
(the name GUA is arbitrary) or use its SMILES code within double quotes
```
python -m auto_martini --smi "N1=C(N)NN=C1N" --mol GUA --top GUA.itp
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
python -m auto_martini --smi "N1=C(N)NN=C1N" --mol GUA --cg gua.gro --top GUA.itp
```
Atomistic coordinates can be written using the `--aa output.gro` option.

## Caveats

For frequently encountered problems, see [FEP](FEP.md).


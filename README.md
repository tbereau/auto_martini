Auto_MARTINI
============
***
*Author:* Tristan Bereau (Max Planck Institute for Polymer Research, Mainz, Germany)  
*Created:* 2014  
***
Automated MARTINI mapping and parametrization of small organic molecules.

## Publication
For a detailed account of the software, see:

Bereau and Kremer, *J Chem Theory Comput*, DOI:10.1021/acs.jctc.5b00056 (2015) [link](http://dx.doi.org/10.1021/acs.jctc.5b00056)

[![DOI for Citing auto_martini](https://img.shields.io/badge/DOI-10.1021%2Facs.jctc.5b00056-blue.svg)](http://dx.doi.org/10.1021/acs.jctc.5b00056)

Please consider citing the paper if you find `auto_martini` useful in your research.

##Installation
`auto-martini` is a python script that requires a number of dependencies:
* `numpy`: see http://docs.scipy.org/doc/numpy/user/install.html
* `rdkit`: see http://www.rdkit.org/docs/Install.html
* `beautifulsoup`: see http://www.crummy.com/software/BeautifulSoup/
* `requests`: see http://docs.python-requests.org/en/latest/user/install/

`numpy` and `rdkit` can be installed by some package managers. Otherwise you'll have to compile it from source. `beautifulsoup` and `requests` can easily be installed using [pip](https://pip.pypa.io/en/latest/) or [easy_install](https://pypi.python.org/pypi/setuptools). In case you do not have root access to your computer to install new software, have a look at [virtualenv](https://pypi.python.org/pypi/virtualenv).

Once all the dependencies are correctly installed, a call to the program should return a usage-information message similar to the following:
```
usage: auto-martini [-h] (--sdf SDF | --smi SMI) --mol MOLNAME [--xyz XYZ]
                    [--gro GRO] [--verbose] [--fpred]
auto-martini: error: argument --mol is required
```

##Usage
To coarse-grain a molecule, simply provide its SMILES code (option `--smi SMI`) or a .SDF file (option `'--sdf file.sdf`). You also need to provide a name for the CG molecule (not longer than 5 characters) using the `--mol` option.  For instance, to coarse grain [guanazole](http://pubchem.ncbi.nlm.nih.gov/summary/summary.cgi?cid=15078), you can either obtain/generate (e.g., from Open Babel) an SDF file:
```
auto-martini --sdf guanazole.sdf --mol GUA
```
(the name GUA is arbitrary) or use its SMILES code within double quotes
```
auto-martini --smi "N1=C(N)NN=C1N" --mol GUA
```
In case no problem arises, it will output the gromacs .itp file:
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
auto-martini --smi "N1=C(N)NN=C1N" --mol GUA --gro gua.gro
```
Atomistic coordinates can be output in XYZ format using the `--xyz output.xyz` option.

##Caveats

###Prediction algorithm

Since ALOGPS, the prediction algorithm for octanol/water partitioning, relies on whole fragments rather than individual atoms, the prediction of certain fragments can pose problem, e.g., small inorganic groups. In this case, `auto-martini` tries to parametrize alternative mappings. If none of them shows successful, the code will return an error.
```
; ERROR: no successful mapping found.
; Try running with the '--fpred' and/or '--verbose' options.
```
As mentioned in the error message, an alternative solution consists of relying on an atom-based partitioning coefficient prediction algorithm (Wildman-Crippen), which is less accurate but can predict any fragment.  In case the `--fpred` option is selected, only fragments for which ALOGPS fail will be predicted using Wildman-Crippen.

###Boost error

Some versions of Boost will fail to correctly exit at the end of the program, generating such output messages:
```
python: /usr/include/boost/thread/pthread/mutex.hpp:108: boost::mutex::~mutex(): Assertion `!posix::pthread_mutex_destroy(&m)' failed.
[1]    31433 abort (core dumped)  ./auto_martini --smi "N1=C(N)NN=C1N" --mol GUA
```
the results provided by the code are unaffected by this error message. Simply ignore it.

###RDKit outdated

Older RDKit versions will report the following error:
```
[...]
distBdAt = Chem.rdMolTransforms.GetBondLength(conf,i,beadId)
AttributeError: 'module' object has no attribute 'GetBondLength'
```
Simply update your version of RDKit. Most package managers will allow you to do this, unless you've installed RDKit from source.

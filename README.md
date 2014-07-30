Auto-MARTINI
============
***
*Author:* Tristan Bereau (Max Planck Institute for Polymer Research, Mainz, Germany)  
*Created:* 2014  
***
Automatized MARTINI mapping and parametrization of small organic molecules.

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

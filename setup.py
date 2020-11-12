'''
Created on March 17, 2019 by Andrew Abi-Mansour

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

import sys, os

# Extract metadata from simulation._version
with open(os.path.join('src', 'auto_martini', '_version.py'), 'r') as fp:
        for line in fp.readlines():
                if '__version__' in line:
                        __version__ = line.split('=')[-1].strip().strip("''")

class Installer(object):

  def __init__(self, repo):
    self.repo = repo
    self.name = repo.split('/')[-1]

  def find(self, fname, path):

    for root, dirs, files in os.walk(path):
      if fname in files:
        return os.path.join(root, fname)

    return None

  def __enter__(self):

    os.system('git clone ' + self.repo)

    os.chdir(self.name)
    os.mkdir('build')
    os.chdir('build')

    python_version = str(sys.version_info[0]) + '.' + str(sys.version_info[1])

    python_lib = self.find('libpython{}.so'.format(python_version), '/') # unix-only 
    python_exec = sys.executable

    if not python_lib:
      print('Could not find any installed python-dev (libpython{}.so).'.format(python_version))
      print('Proceeding ...')
      cm_args = ' -DPYTHON_EXECUTABLE={} -DCMAKE_INSTALL_PREFIX=$HOME/.local'.format(python_exec)
    else:
      cm_args = ' -DPYTHON_LIBRARY={} -DPYTHON_EXECUTABLE={} -DCMAKE_INSTALL_PREFIX=$HOME/.local'.format(python_lib, python_exec)

    os.system('cmake .. ' + cm_args)

    if 'RDKITPROC' in os.environ:
      cmd = 'make install -j{}'.format(os.environ['RDKITPROC'])
    else:
      import multiprocessing
      cmd = 'make install -j{}'.format(multiprocessing.cpu_count())

    os.system(cmd)

  def __exit__(self, *a):
    os.chdir(os.path.join('..','..'))
    sys.path.append(os.path.join(os.getcwd(), self.name, 'lib'))

# check if rdkit is installed ... else compile it from source
try:
  import rdkit
except Exception:
  print('rdkit not found. Attempting to compile rdkit from source ...')
  with Installer(repo='https://github.com/rdkit/rdkit') as _:
    pass

if __name__ == '__main__':

  from setuptools import setup, find_packages

  try:
  	from Cython.Build import cythonize
  	import numpy
  	optimal_list = cythonize("src/auto_martini/optimization.pyx")
  	include_dirs = [numpy.get_include()]
  except:
    print('Failed to cythonize optimization module. For optimal performance, make sure Cython is properly installed.')
    optimal_list = []
    include_dirs = []

  setup(
      name = "auto_martini",
      version = __version__,
      author = ["Tristan Bereau", "Andrew Abi-Mansour"],
      author_email = ["bereau [at] mpip-mainz.mpg.de", "andrew.gaam [at] gmail.com"],
      description = ("A tool for automatic MARTINI mapping and parametrization of small organic molecules "),
      license = "GPL v2",
      keywords = "Coarse-grained Molecular Dynamics, MARTINI Force Field",
      url = "https://github.com/Andrew-AbiMansour/Auto_MARTINI",
      packages=find_packages('src'),
      package_dir={'auto_martini':'src/auto_martini'},
      package_data={'test': ['test/*.sdf'],},
      include_package_data=True,
      install_requires=['numpy', 'bs4', 'pytool', 'lxml', 'requests'],
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

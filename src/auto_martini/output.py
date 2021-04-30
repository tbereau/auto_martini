'''
Created on March 13, 2019 by Andrew Abi-Mansour

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

from .common import *
logger = logging.getLogger(__name__)

def output_gro(sites, site_names, molname):
    """Output GRO file of CG structure"""
    logger.debug('Entering output_gro()')
    num_beads = len(sites)
    gro_out = ""
    if len(sites) != len(site_names):
        logger.warning('Error. Incompatible number of beads and bead names.')
        exit(1)
    gro_out += "{:s} generated from auto_martini\n".format(molname)
    gro_out += "{:5d}\n".format(num_beads)
    for i in range(num_beads):
        gro_out += "{:5d}{:<6s} {:3s}{:5d}{:8.3f}{:8.3f}{:8.3f}\n".format(
            i + 1, molname, site_names[i], i + 1, sites[i][0] / 10.,
            sites[i][1] / 10., sites[i][2] / 10.)
    gro_out += "{:10.5f}{:10.5f}{:10.5f}\n".format(10., 10., 10.)
    return gro_out

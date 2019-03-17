# !/usr/bin/python
#  -*- coding: utf8 -*-

'''
  Created on March 13, 2019 by Andrew Abi-Mansour

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
'''

from auto_martini.common import *

def output_gro(output_file, sites, site_names, molname):
    """Output GRO file of CG structure"""
    logging.debug('Entering output_gro()')
    num_beads = len(sites)

    if len(sites) != len(site_names):
        logging.warning('Error. Incompatible number of beads and bead names.')
        exit(1)
    if output_file[-4:] != ".gro":
        output_file += ".gro"
    try:
        with open(output_file, 'w') as f:
            f.write("{:s} generated from {:s}\n".format(
                molname, os.path.basename(__file__)))
            f.write("{:5d}\n".format(num_beads))
            for i in range(num_beads):
                f.write("{:5d}{:<6s} {:3s}{:5d}{:8.3f}{:8.3f}{:8.3f}\n".format(
                    i + 1, molname, site_names[i], i + 1, sites[i][0] / 10.,
                    sites[i][1] / 10., sites[i][2] / 10.))
            f.write("{:10.5f}{:10.5f}{:10.5f}\n".format(10., 10., 10.))
            f.close()
    except IOError:
        logging.warning('Can\'t write to file %s' % output_file)
        exit(1)
    return
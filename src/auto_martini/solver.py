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

  BUGS found and fixed:

  - function substruct2smi() returned SMILES string in lower case letters
  - function printAtoms() uses undefined objects "hbondA" and "hbondD"
  - function checkAdditivity() defined "mad" variable which is already defined as a function
  - function genMoleculeSDF() did not use "sdf" input argument
  - function printBonds() used undefined "atomPartitioning" object
  - function Substructure() no longer creates then reads "tmp-auto-martini.smi" file


  TODO: make this run in Python 3

'''

from auto_martini.common import *

from auto_martini import output
from auto_martini import topology
from auto_martini import optimization

def get_coords(conformer, sites, avg_pos, ringatoms_flat):
    """Extract coordinates of CG beads"""
    logging.debug('Entering get_coords()')
    # CG beads are averaged over best trial combinations for all
    # non-aromatic atoms.
    logging.debug('Entering get_coords()')
    site_coords = []
    for i in range(len(sites)):
        if sites[i] in ringatoms_flat:
            site_coords.append(np.array([conformer.GetAtomPosition(int(sites[i]))[j] for j in range(3)]))
        else:
            # Use average
            site_coords.append(np.array(avg_pos[i]))
    return site_coords


def check_additivity(forcepred, beadtypes, molecule):
    """Check additivity assumption between sum of free energies of CG beads
    and free energy of whole molecule"""
    logging.debug('Entering check_additivity()')
    # If there's only one bead, don't check.
    if len(beadtypes) == 1:
        return True
    sum_frag = 0.0
    rings = False
    logging.info('; Bead types: %s' % beadtypes)
    for bead in beadtypes:
        if bead[0] == "S":
            bead = bead[1:]
            rings = True
        delta_f_types = topology.read_delta_f_types()
        sum_frag += delta_f_types[bead]
    # Wildman-Crippen log_p
    wc_log_p = rdMolDescriptors.CalcCrippenDescriptors(molecule)[0]
    # Write out SMILES string of entire molecule
    tmpfile = 'tmp-auto-martini.smi'
    sw = rdmolfiles.SmilesWriter(tmpfile, nameHeader='')
    rdmolfiles.SmilesWriter.write(sw, molecule)
    sw.close()
    # Read file
    s = ''
    try:
        f = open(tmpfile, 'r')
        s = f.readlines()
        f.close()
    except IOError as e:
        logging.warning('Error. Can\'t read file %s' % tmpfile)
        logging.warning(e)
        exit(1)
    os.remove(tmpfile)
    whole_mol_dg = topology.smi2alogps(forcepred, s[1].split()[0], wc_log_p, "MOL", True)
    m_ad = math.fabs((whole_mol_dg - sum_frag) / whole_mol_dg)
    logging.info('; Mapping additivity assumption ratio: %7.4f (whole vs sum: %7.4f vs. %7.4f)'
                % (m_ad, whole_mol_dg, sum_frag))
    if (not rings and m_ad < 0.5) or rings:
        return True
    else:
        return False

def cg_molecule(molecule, molname, aa_output=None, cg_output=None, forcepred=False):
    """Main routine to coarse-grain molecule"""
    # Get molecule's features

    feats = topology.extract_features(molecule)

    # Get list of heavy atoms and their coordinates
    list_heavy_atoms, list_heavyatom_names = topology.get_atoms(molecule)
    conf, heavy_atom_coords = topology.get_heavy_atom_coords(molecule)

    # Identify ring-type atoms
    ring_atoms = topology.get_ring_atoms(molecule)

    # Get Hbond information
    hbond_a = topology.get_hbond_a(feats)
    hbond_d = topology.get_hbond_d(feats)

    # Flatten list of ring atoms
    ring_atoms_flat = list(chain.from_iterable(ring_atoms))

    # Optimize coarse-grained bead positions -- keep all possibilities in case something goes
    # wrong later in the code.
    list_cg_beads, list_bead_pos = optimization.find_bead_pos(molecule, conf, list_heavy_atoms, heavy_atom_coords, ring_atoms,
                                                 ring_atoms_flat)

    # Loop through best 1% cg_beads and avg_pos
    cg_bead_names = []
    cg_bead_coords = []
    max_attempts = int(math.ceil(0.5 * len(list_cg_beads)))
    logging.info('; Max. number of attempts: %s' % max_attempts)
    attempt = 0

    while attempt < max_attempts:
        cg_beads = list_cg_beads[attempt]
        bead_pos = list_bead_pos[attempt]
        success = True

        # Extract position of coarse-grained beads
        cg_bead_coords = get_coords(conf, cg_beads, bead_pos, ring_atoms_flat)

        # Partition atoms into coarse-grained beads
        atom_partitioning = optimization.voronoi_atoms(cg_bead_coords, heavy_atom_coords)
        logging.info('; Atom partitioning: %s' % atom_partitioning)

        cg_bead_names, bead_types = topology.print_atoms(molname, forcepred, cg_beads, molecule, hbond_a, hbond_d, atom_partitioning, ring_atoms, ring_atoms_flat, True)

        if not cg_bead_names:
            success = False
        # Check additivity between fragments and entire molecule
        if not check_additivity(forcepred, bead_types, molecule):
            success = False
        # Bond list
        try:
            bond_list, const_list = topology.print_bonds(cg_beads, molecule, atom_partitioning, cg_bead_coords, ring_atoms, True)
        except (NameError, ValueError):
            success = False

        if success:
            topology.print_header(molname)
            cg_bead_names, bead_types = topology.print_atoms(molname, forcepred, cg_beads, molecule, hbond_a, hbond_d,
                                                    atom_partitioning, ring_atoms, ring_atoms_flat, False)

            bond_list, const_list = topology.print_bonds(cg_beads, molecule, atom_partitioning, cg_bead_coords, ring_atoms,
                                                False)
            topology.print_angles(cg_beads, molecule, atom_partitioning, cg_bead_coords, bond_list, const_list, ring_atoms)
            topology.print_dihedrals(cg_beads, const_list, ring_atoms, cg_bead_coords)

            # We've reached all the way here, exit the while loop
            attempt = max_attempts + 1
        else:
            attempt += 1
            logging.info('------------------------------------------------------')

    if attempt == max_attempts:
        err = "; ERROR: no successful mapping found.\n" + \
              "; Try running with the '--fpred' and/or '--verbose' options.\n"
        sys.stderr.write(err)
        exit(1)

    # Optional all-atom output to GRO file
    if aa_output:
        output.output_gro(aa_output, heavy_atom_coords, list_heavyatom_names, "MOL")

    # Optional coarse-grained output to GRO file
    if cg_output:
        output.output_gro(cg_output, cg_bead_coords, cg_bead_names, molname)

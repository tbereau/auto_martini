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
from .common import __version__

from . import output
from . import topology
from . import optimization

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

    for bead in beadtypes:
        if bead[0] == "S":
            bead = bead[1:]
            rings = True
        delta_f_types = topology.read_delta_f_types()
        sum_frag += delta_f_types[bead]
    # Wildman-Crippen log_p
    wc_log_p = rdMolDescriptors.CalcCrippenDescriptors(molecule)[0]
    
    whole_mol_dg = topology.smi2alogps(forcepred, Chem.MolToSmiles(molecule), wc_log_p, "MOL", True)
    m_ad = math.fabs((whole_mol_dg - sum_frag) / whole_mol_dg)
    logging.info('Mapping additivity assumption ratio: %7.4f (whole vs sum: %7.4f vs. %7.4f)'
                % (m_ad, whole_mol_dg, sum_frag))

    if (not rings and m_ad < 0.5) or rings:
        return True
    else:
        return False

def cg_molecule(molecule, molname, topfname, aa_output=None, cg_output=None, forcepred=False):
    """Main routine to coarse-grain molecule"""
    # Get molecule's features

    logging.info('Running Auto_Martini v{}'.format(__version__))

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
    logging.info('Max. number of attempts: %s' % max_attempts)
    attempt = 1

    while attempt < max_attempts:
        cg_beads = list_cg_beads[attempt]
        bead_pos = list_bead_pos[attempt]
        success = True

        # Extract position of coarse-grained beads
        cg_bead_coords = get_coords(conf, cg_beads, bead_pos, ring_atoms_flat)

        # Partition atoms into coarse-grained beads
        atom_partitioning = optimization.voronoi_atoms(cg_bead_coords, heavy_atom_coords)
        logging.info('; Atom partitioning: %s' % atom_partitioning)

        cg_bead_names, bead_types, _ = topology.print_atoms(molname, forcepred, cg_beads, molecule, hbond_a, hbond_d, atom_partitioning, ring_atoms, ring_atoms_flat, True)

        if not cg_bead_names:
            success = False
        # Check additivity between fragments and entire molecule
        if not check_additivity(forcepred, bead_types, molecule):
            success = False
        # Bond list
        try:
            bond_list, const_list, _ = topology.print_bonds(cg_beads, molecule, atom_partitioning, cg_bead_coords, ring_atoms, trial=True)
        except:
            raise

        if success:
            header_write = topology.print_header(molname)
            cg_bead_names, bead_types, atoms_write = topology.print_atoms(molname, forcepred, cg_beads, molecule, hbond_a, hbond_d,
                                                    atom_partitioning, ring_atoms, ring_atoms_flat, trial=False)

            bond_list, const_list, bonds_write = topology.print_bonds(cg_beads, molecule, atom_partitioning, cg_bead_coords, ring_atoms,
                                                False)
            angles_write = topology.print_angles(cg_beads, molecule, atom_partitioning, cg_bead_coords, bond_list, const_list, ring_atoms)
            dihedrals_write = topology.print_dihedrals(cg_beads, const_list, ring_atoms, cg_bead_coords)

            with open(topfname, 'w') as fp:
                fp.write(header_write + atoms_write + bonds_write +  angles_write + dihedrals_write) 

            print('Converged to solution in {} iteration(s)'.format(attempt))
            break
        else:
            attempt += 1

    if attempt == max_attempts:
        raise RuntimeError('ERROR: no successful mapping found.\nTry running with the --fpred and/or --verbose options.')

    # Optional all-atom output to GRO file
    if aa_output:
        output.output_gro(aa_output, heavy_atom_coords, list_heavyatom_names, "MOL")

    # Optional coarse-grained output to GRO file
    if cg_output:
        output.output_gro(cg_output, cg_bead_coords, cg_bead_names, molname)

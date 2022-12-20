"""
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
"""

from . import optimization, output, topology
from .common import *

logger = logging.getLogger(__name__)


def get_coords(conformer, sites, avg_pos, ringatoms_flat):
    """Extract coordinates of CG beads"""
    logger.debug("Entering get_coords()")
    # CG beads are averaged over best trial combinations for all
    # non-aromatic atoms.
    logger.debug("Entering get_coords()")
    site_coords = []
    for i in range(len(sites)):
        if sites[i] in ringatoms_flat:
            site_coords.append(
                np.array([conformer.GetAtomPosition(int(sites[i]))[j] for j in range(3)])
            )
        else:
            # Use average
            site_coords.append(np.array(avg_pos[i]))
    return site_coords


def check_additivity(forcepred, beadtypes, molecule):
    """Check additivity assumption between sum of free energies of CG beads
    and free energy of whole molecule"""
    logger.debug("Entering check_additivity()")
    # If there's only one bead, don't check.
    sum_frag = 0.0
    rings = False
    logger.info("; Bead types: %s" % beadtypes)
    for bead in beadtypes:
        if bead[0] == "S":
            bead = bead[1:]
            rings = True
        delta_f_types = topology.read_delta_f_types()
        sum_frag += delta_f_types[bead]
    # Wildman-Crippen log_p
    wc_log_p = rdMolDescriptors.CalcCrippenDescriptors(molecule)[0]
    # Get SMILES string of entire molecule

    s = Chem.MolToSmiles(molecule)

    whole_mol_dg = topology.smi2alogps(forcepred, s, wc_log_p, "MOL", True)
    if whole_mol_dg != 0:
        m_ad = math.fabs((whole_mol_dg - sum_frag) / whole_mol_dg)
        logger.info(
            "; Mapping additivity assumption ratio: %7.4f (whole vs sum: %7.4f vs. %7.4f)"
            % (m_ad, whole_mol_dg / (-4.184), sum_frag / (-4.184))
        )
        if len(beadtypes) == 1:
            return True
        if (not rings and m_ad < 0.5) or rings:
            return True
        else:
            return False
    else:
        return False


class Cg_molecule:
    """Main class to coarse-grain molecule"""

    def __init__(self, molecule, molname, topfname=None, forcepred=False):
        self.heavy_atom_coords = None
        self.list_heavyatom_names = None
        self.atom_partitioning = None
        self.cg_bead_names = []
        self.cg_bead_coords = []
        self.topout = None

        # Get molecule's features

        logger.info("Entering cg_molecule()")

        feats = topology.extract_features(molecule)

        # Get list of heavy atoms and their coordinates
        list_heavy_atoms, self.list_heavyatom_names = topology.get_atoms(molecule)

        conf, self.heavy_atom_coords = topology.get_heavy_atom_coords(molecule)

        # Identify ring-type atoms
        ring_atoms = topology.get_ring_atoms(molecule)

        # Get Hbond information
        hbond_a = topology.get_hbond_a(feats)
        hbond_d = topology.get_hbond_d(feats)

        # Flatten list of ring atoms
        ring_atoms_flat = list(chain.from_iterable(ring_atoms))

        # Optimize coarse-grained bead positions -- keep all possibilities in case something goes
        # wrong later in the code.
        list_cg_beads, list_bead_pos = optimization.find_bead_pos(
            molecule,
            conf,
            list_heavy_atoms,
            self.heavy_atom_coords,
            ring_atoms,
            ring_atoms_flat,
        )

        # Loop through best 1% cg_beads and avg_pos
        max_attempts = int(math.ceil(0.5 * len(list_cg_beads)))
        logger.info(f"Max. number of attempts: {max_attempts}")
        attempt = 0

        while attempt < max_attempts:
            cg_beads = list_cg_beads[attempt]
            bead_pos = list_bead_pos[attempt]
            success = True

            # Remove mappings with bead numbers less than most optimal mapping.
            if (
                len(cg_beads) < len(list_cg_beads[0])
                and (len(list_heavy_atoms) - (5 * len(cg_beads))) > 3
            ):
                success = False

            # Extract position of coarse-grained beads
            self.cg_bead_coords = get_coords(conf, cg_beads, bead_pos, ring_atoms_flat)

            # Partition atoms into coarse-grained beads
            self.atom_partitioning = optimization.voronoi_atoms(
                self.cg_bead_coords, self.heavy_atom_coords
            )
            logger.info("; Atom partitioning: {atom_partitioning}")

            self.cg_bead_names, bead_types, _ = topology.print_atoms(
                molname,
                forcepred,
                cg_beads,
                molecule,
                hbond_a,
                hbond_d,
                self.atom_partitioning,
                ring_atoms,
                ring_atoms_flat,
                True,
            )

            if not self.cg_bead_names:
                success = False
            # Check additivity between fragments and entire molecule
            if not check_additivity(forcepred, bead_types, molecule):
                success = False
            # Bond list
            try:
                bond_list, const_list, _ = topology.print_bonds(
                    cg_beads,
                    molecule,
                    self.atom_partitioning,
                    self.cg_bead_coords,
                    ring_atoms,
                    trial=True,
                )
            except Exception:
                raise

            # I added errval below from the master branch ... not sure where to use this anywhere, possibly leave for debugging
            if not ring_atoms and (len(bond_list) + len(const_list)) >= len(self.cg_bead_names):
                errval = 3
                success = False
            if (len(bond_list) + len(const_list)) < len(self.cg_bead_names) - 1:
                errval = 5
                success = False
            if len(cg_beads) != len(self.cg_bead_names):
                success = False
                errval = 8

            if success:
                header_write = topology.print_header(molname)
                self.cg_bead_names, bead_types, atoms_write = topology.print_atoms(
                    molname,
                    forcepred,
                    cg_beads,
                    molecule,
                    hbond_a,
                    hbond_d,
                    self.atom_partitioning,
                    ring_atoms,
                    ring_atoms_flat,
                    trial=False,
                )

                bond_list, const_list, bonds_write = topology.print_bonds(
                    cg_beads,
                    molecule,
                    self.atom_partitioning,
                    self.cg_bead_coords,
                    ring_atoms,
                    False,
                )
                angles_write, angle_list = topology.print_angles(
                    cg_beads,
                    molecule,
                    self.atom_partitioning,
                    self.cg_bead_coords,
                    bond_list,
                    const_list,
                    ring_atoms,
                )

                if not angles_write and len(bond_list) > 1:
                    errval = 2
                if bond_list and angle_list:
                    if (len(bond_list) + len(const_list)) < 2 and len(angle_list) > 0:
                        errval = 6
                    if (
                        not ring_atoms
                        and (len(bond_list) + len(const_list)) - len(angle_list) != 1
                    ):
                        errval = 7

                dihedrals_write = topology.print_dihedrals(
                    cg_beads, const_list, ring_atoms, self.cg_bead_coords
                )

                self.topout = (
                    header_write + atoms_write + bonds_write + angles_write + dihedrals_write
                )
                if topfname:
                    with open(topfname, "w") as fp:
                        fp.write(self.topout)

                print("Converged to solution in {} iteration(s)".format(attempt + 1))
                break
            else:
                attempt += 1

        if attempt == max_attempts:
            raise RuntimeError(
                "ERROR: no successful mapping found.\nTry running with the --fpred and/or --verbose options."
            )

    def output_aa(self, aa_output=None, molname="MOL"):
        # Optional all-atom output to GRO file
        aa_out = output.output_gro(self.heavy_atom_coords, self.list_heavyatom_names, molname)
        if aa_output:
            with open(aa_output, "w") as fp:
                fp.write(aa_out)
        else:
            return aa_out

    def output_cg(self, cg_output=None, molname="MOL"):
        # Optional coarse-grained output to GRO file
        cg_out = output.output_gro(self.cg_bead_coords, self.cg_bead_names, molname)
        if cg_output:
            with open(cg_output, "w") as fp:
                fp.write(cg_out)
        else:
            return cg_out

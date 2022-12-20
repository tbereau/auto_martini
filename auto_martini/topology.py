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

from sys import exit

from auto_martini._version import __version__

from .common import *

logger = logging.getLogger(__name__)

# For feature extraction
fdefName = os.path.join(RDConfig.RDDataDir, "BaseFeatures.fdef")
factory = ChemicalFeatures.BuildFeatureFactory(fdefName)


def read_delta_f_types():
    """Returns delta_f types dictionary
    Measured octanol/water free energies from MARTINI
    Data used TI, not partitioning of Marrink et al. JPCB 2007"""
    delta_f_types = dict()
    delta_f_types["Qda"] = -15.04
    delta_f_types["Qa"] = -15.04
    delta_f_types["Qd"] = -15.04
    delta_f_types["Q0"] = -22.35
    delta_f_types["P5"] = -8.88
    delta_f_types["P4"] = -9.30
    delta_f_types["P3"] = -8.81
    delta_f_types["P2"] = -3.85
    delta_f_types["P1"] = -2.26
    delta_f_types["Nda"] = 2.49
    delta_f_types["Na"] = 2.49
    delta_f_types["Nd"] = 2.49
    delta_f_types["N0"] = 4.22
    delta_f_types["C5"] = 6.93
    delta_f_types["C4"] = 10.14
    delta_f_types["C3"] = 12.26
    delta_f_types["C2"] = 13.74
    delta_f_types["C1"] = 14.20
    return delta_f_types


def gen_molecule_smi(smi):
    """Generate mol object from smiles string"""
    logger.debug("Entering gen_molecule_smi()")
    errval = 0
    if "." in smi:
        logger.warning("Error. Only one molecule may be provided.")
        logger.warning(smi)
        errval = 4
        exit(1)
    # If necessary, adjust smiles for Aromatic Ns
    # Redirect current stderr in log file
    stderr_fd = None
    stderr_save = None
    try:
        stderr_fileno = sys.stderr.fileno()
        stderr_save = os.dup(stderr_fileno)
        stderr_fd = open("sanitize.log", "w")
        os.dup2(stderr_fd.fileno(), stderr_fileno)
    except Exception:
        stderr_fileno = None
    # Get smiles without sanitization
    molecule = Chem.MolFromSmiles(smi, False)
    try:
        cp = Chem.Mol(molecule)
        Chem.SanitizeMol(cp)
        # Close log file and restore old sys err
        if stderr_fileno is not None:
            stderr_fd.close()
            os.dup2(stderr_save, stderr_fileno)
        molecule = cp
    except ValueError:
        logger.warning("Bad smiles format %s found" % smi)
        nm = AdjustAromaticNs(molecule)
        if nm is not None:
            Chem.SanitizeMol(nm)
            molecule = nm
            smi = Chem.MolToSmiles(nm)
            logger.warning("Fixed smiles format to %s" % smi)
        else:
            logger.warning("Smiles cannot be adjusted %s" % smi)
            errval = 1
    # Continue
    molecule = Chem.AddHs(molecule)
    AllChem.EmbedMolecule(
        molecule, randomSeed=1, useRandomCoords=True
    )  # Set Seed for random coordinate generation = 1.
    try:
        AllChem.UFFOptimizeMolecule(molecule)
    except ValueError as e:
        logger.warning("%s" % e)
        exit(1)
    return molecule, errval


def gen_molecule_sdf(sdf):
    """Generate mol object from SD file"""
    logger.debug("Entering gen_molecule_sdf()")
    suppl = Chem.SDMolSupplier(sdf)
    if len(suppl) > 1:
        print("Error. Only one molecule may be provided.")
        exit(1)
    molecule = ""
    for molecule in suppl:
        if molecule is None:
            print("Error. Can't read molecule.")
            exit(1)
    return molecule


def print_header(molname):
    """Print topology header"""
    logger.debug("Entering print_header()")

    text = "; GENERATED WITH auto_Martini v{} for {}\n".format(__version__, molname)

    info = (
        "; Developed by: Kiran Kanekal, Tristan Bereau, and Andrew Abi-Mansour\n\n"
        + "[moleculetype]\n"
        + "; molname       nrexcl\n"
        + "  {:5s}         2\n\n".format(molname)
        + "[atoms]\n"
        + "; id    type    resnr   residue  atom    cgnr    charge  smiles\n"
    )

    return text + info


def letter_occurrences(string):
    """Count letter occurences"""
    logger.debug("Entering letter_occurrences()")
    frequencies = defaultdict(lambda: 0)
    for character in string:
        if character.isalnum():
            frequencies[character.upper()] += 1
    return frequencies


def get_charge(molecule):
    """Get net charge of molecule"""
    logger.debug("Entering get_charge()")
    return Chem.rdmolops.GetFormalCharge(molecule)


def get_hbond_a(features):
    """Get Hbond acceptor information"""
    logger.debug("Entering get_hbond_a()")
    hbond = []
    for feat in features:
        if feat.GetFamily() == "Acceptor":
            for i in feat.GetAtomIds():
                if i not in hbond:
                    hbond.append(i)
    return hbond


def get_hbond_d(features):
    """Get Hbond donor information"""
    logger.debug("Entering get_hbond_d()")
    hbond = []
    for feat in features:
        if feat.GetFamily() == "Donor":
            for i in feat.GetAtomIds():
                if i not in hbond:
                    hbond.append(i)
    return hbond


def get_atoms(molecule):
    """List all heavy atoms"""
    logger.debug("Entering get_atoms()")
    conformer = molecule.GetConformer()
    num_atoms = conformer.GetNumAtoms()
    list_heavyatoms = []
    list_heavyatomnames = []

    atoms = np.arange(num_atoms)
    for i in np.nditer(atoms):
        atom_name = molecule.GetAtomWithIdx(int(atoms[i])).GetSymbol()
        if atom_name != "H":
            list_heavyatoms.append(atoms[i])
            list_heavyatomnames.append(atom_name)

    if len(list_heavyatoms) == 0:
        print("Error. No heavy atom found.")
        exit(1)
    return list_heavyatoms, list_heavyatomnames


def get_ring_atoms(molecule):
    """Get ring atoms"""
    logger.debug("Entering get_ring_atoms()")
    ringatoms = []
    ringinfo = molecule.GetRingInfo()
    rings = ringinfo.AtomRings()
    for at in rings:
        ring = list(sorted(at))
        ringatoms.append(ring)
    return ringatoms


def get_heavy_atom_coords(molecule):
    """Extract atomic coordinates of heavy atoms in molecule mol"""
    logger.debug("Entering get_heavy_atom_coords()")
    heavyatom_coords = []
    conformer = molecule.GetConformer()
    # number of atoms in mol
    num_atoms = molecule.GetConformer().GetNumAtoms()
    for i in range(num_atoms):
        if molecule.GetAtomWithIdx(i).GetSymbol() != "H":
            heavyatom_coords.append(np.array([conformer.GetAtomPosition(i)[j] for j in range(3)]))

    return conformer, heavyatom_coords


def extract_features(molecule):
    """Extract features of molecule"""
    logger.debug("Entering extract_features()")
    features = factory.GetFeaturesForMol(molecule)
    return features


def substruct2smi(molecule, partitioning, cg_bead, cgbeads, ringatoms):
    """Substructure to smiles conversion; also output Wildman-Crippen log_p;
    and charge of group."""

    frag = rdchem.EditableMol(molecule)

    # fragment smi: [H]N([H])c1nc(N([H])[H])n([H])n1

    num_atoms = molecule.GetConformer().GetNumAtoms()
    # First delete all hydrogens
    for i in range(num_atoms):
        if molecule.GetAtomWithIdx(i).GetSymbol() == "H":
            # find atom from coordinates
            submol = frag.GetMol()
            for j in range(submol.GetConformer().GetNumAtoms()):
                if (
                    molecule.GetConformer().GetAtomPosition(i)[0]
                    == submol.GetConformer().GetAtomPosition(j)[0]
                ):
                    frag.RemoveAtom(j)
    # Identify atoms involved in same ring as cg_bead (only one ring)
    atoms_in_ring = []
    for ring in ringatoms:
        if cgbeads[cg_bead] in ring:
            atoms_in_ring = ring[:]  # CHANGED
            break
    # Then heavy atoms that aren't part of the CG bead (except those
    # involved in the same ring).
    for i in partitioning.keys():
        if partitioning[i] != cg_bead and i not in atoms_in_ring:
            # find atom from coordinates
            submol = frag.GetMol()
            for j in range(submol.GetConformer().GetNumAtoms()):
                if (
                    molecule.GetConformer().GetAtomPosition(i)[0]
                    == submol.GetConformer().GetAtomPosition(j)[0]
                ):
                    frag.RemoveAtom(j)
    # Wildman-Crippen log_p
    wc_log_p = rdMolDescriptors.CalcCrippenDescriptors(frag.GetMol())[0]

    # Charge -- look at atoms that are only part of the bead (no ring rule)
    chg = 0
    for i in partitioning.keys():
        if partitioning[i] == cg_bead:
            chg += molecule.GetAtomWithIdx(i).GetFormalCharge()

    smi = Chem.MolToSmiles(Chem.rdmolops.AddHs(frag.GetMol(), addCoords=True))

    # fragment smi: Nc1ncnn1 ---------> FAILURE! Need to fix this Andrew! For now, just a hackish soln:
    # smi = smi.lower() if smi.islower() else smi.upper()

    return smi, wc_log_p, chg


def print_atoms(
    molname,
    forcepred,
    cgbeads,
    molecule,
    hbonda,
    hbondd,
    partitioning,
    ringatoms,
    ringatoms_flat,
    trial=False,
):
    """Print CG Atoms in itp format"""

    logger.debug("Entering print_atoms()")
    atomnames = []
    beadtypes = []
    text = ""

    for bead in range(len(cgbeads)):
        # Determine SMI of substructure
        try:
            smi_frag, wc_log_p, charge = substruct2smi(
                molecule, partitioning, bead, cgbeads, ringatoms
            )
        except Exception:
            raise

        atom_name = ""
        for character, count in sorted(six.iteritems(letter_occurrences(smi_frag))):
            try:
                float(character)
            except ValueError:
                if count == 1:
                    atom_name += "{:s}".format(character)
                else:
                    atom_name += "{:s}{:s}".format(character, str(count))

        # Get charge for smi_frag
        mol_frag, errval = gen_molecule_smi(smi_frag)
        charge_frag = get_charge(mol_frag)

        if errval == 0:
            # frag_heavyatom_coord = get_heavy_atom_coords(mol_frag)
            # frag_HA_coord_towrite = frag_heavyatom_coord[1]
            # frag_HA_coord_towrite[:0] = [smi_frag]
            # frag_heavyatom_coord_list.append(frag_HA_coord_towrite)
            charge_frag = get_charge(mol_frag)

            # Extract ALOGPS free energy
            try:
                if charge_frag == 0:
                    alogps = smi2alogps(forcepred, smi_frag, wc_log_p, bead + 1, trial)
                else:
                    alogps = 0.0
            except (NameError, TypeError, ValueError):
                return atomnames, beadtypes, errval

            hbond_a_flag = 0
            for at in hbonda:
                if partitioning[at] == bead:
                    hbond_a_flag = 1
                    break
            hbond_d_flag = 0
            for at in hbondd:
                if partitioning[at] == bead:
                    hbond_d_flag = 1
                    break

            in_ring = cgbeads[bead] in ringatoms_flat
            bead_type = determine_bead_type(alogps, charge, hbond_a_flag, hbond_d_flag, in_ring)
            atom_name = ""
            name_index = 0
            while atom_name in atomnames or name_index == 0:
                name_index += 1
                atom_name = "{:1s}{:02d}".format(bead_type[0], name_index)
            atomnames.append(atom_name)

            if not trial:
                text = (
                    text
                    + "    {:<5d} {:5s}   1     {:5s}     {:7s} {:<5d}    {:2d}   ; {:s}\n".format(
                        bead + 1,
                        bead_type,
                        molname,
                        atom_name,
                        bead + 1,
                        charge,
                        smi_frag,
                    )
                )
            beadtypes.append(bead_type)

    text = text + "\n"

    return atomnames, beadtypes, text


def print_bonds(cgbeads, molecule, partitioning, cgbead_coords, ringatoms, trial=False):
    """print CG bonds in itp format"""
    logger.debug("Entering print_bonds()")

    # Bond information
    bondlist = []
    constlist = []
    text = ""

    if len(cgbeads) > 1:
        for i in range(len(cgbeads)):
            for j in range(i + 1, len(cgbeads)):
                dist = np.linalg.norm(cgbead_coords[i] - cgbead_coords[j]) * 0.1
                if dist < 0.65:
                    # Are atoms part of the same ring
                    in_ring = False
                    for ring in ringatoms:
                        if cgbeads[i] in ring and cgbeads[j] in ring:
                            in_ring = True
                            break
                    if in_ring:
                        constlist.append([i, j, dist])
                    else:
                        # Check that the bond is not too short
                        if dist < 0.2:
                            raise NameError("Bond too short")
                        # Look for a bond between an atom of i and an atom of j
                        found_connection = False
                        atoms_in_bead_i = []
                        for aa in partitioning.keys():
                            if partitioning[aa] == i:
                                atoms_in_bead_i.append(aa)
                        atoms_in_bead_j = []
                        for aa in partitioning.keys():
                            if partitioning[aa] == j:
                                atoms_in_bead_j.append(aa)
                        for ib in range(len(molecule.GetBonds())):
                            abond = molecule.GetBondWithIdx(ib)
                            if (
                                abond.GetBeginAtomIdx() in atoms_in_bead_i
                                and abond.GetEndAtomIdx() in atoms_in_bead_j
                            ) or (
                                abond.GetBeginAtomIdx() in atoms_in_bead_j
                                and abond.GetEndAtomIdx() in atoms_in_bead_i
                            ):
                                found_connection = True
                        if found_connection:
                            bondlist.append([i, j, dist])

        for ring in ringatoms:
            # Only keep one bond between a ring and a given external bead
            for i in range(len(cgbeads)):
                at = cgbeads[i]
                if at not in ring:
                    bonds_to_ring = []
                    for b in bondlist:
                        bead = i
                        if (cgbeads[b[0]] in ring and b[1] == bead) or (
                            b[0] == bead and cgbeads[b[1]] in ring
                        ):
                            bonds_to_ring.append(b)
                    # keep closest
                    closest_bond = [-1, -1, 1000.0]
                    for r in range(len(bonds_to_ring)):
                        if bonds_to_ring[r][2] < closest_bond[2]:
                            closest_bond = bonds_to_ring[r]
                    # Delete the other bonds
                    for b in bonds_to_ring:
                        if b != closest_bond:
                            bondlist.remove(b)
            bead_bonded_to_ring = []
            for i in range(len(cgbeads)):
                atoms_in_bead = []
                for key, val in six.iteritems(partitioning):
                    if val == i:
                        atoms_in_bead.append(key)
                for b in bondlist:
                    if (cgbeads[b[0]] in ring and b[1] in atoms_in_bead) or (
                        b[0] in atoms_in_bead and cgbeads[b[1]] in ring
                    ):
                        bead_bonded_to_ring.append(i)
            # Delete bond between 2 beads if they're both linked
            # to the same ring.
            for i in range(len(cgbeads)):
                for j in range(i + 1, len(cgbeads)):
                    if cgbeads[i] in bead_bonded_to_ring and cgbeads[j] in bead_bonded_to_ring:
                        for b in bondlist:
                            if (b[0] == i and b[1] == j) or (b[0] == j and b[1] == i):
                                bondlist.remove(b)
        # Replace bond by constraint if both atoms have constraints
        # to the same third atom

        bond_list_idx = 0
        while bond_list_idx < len(bondlist):
            b = bondlist[bond_list_idx]
            for i in range(len(cgbeads)):
                const_i = False
                const_j = False
                for c in constlist:
                    if (c[0] == b[0] and c[1] == i) or (c[0] == i and c[1] == b[0]):
                        const_i = True
                    if (c[0] == b[1] and c[1] == i) or (c[0] == i and c[1] == b[1]):
                        const_j = True
                if const_i and const_j:
                    constlist.append(b)
                    bondlist.remove(b)
                    # Start over
                    bond_list_idx = -1
            bond_list_idx += 1

        # Go through list of constraints. If we find an extra
        # possible constraint between beads that have constraints,
        # add it.
        beads_with_const = []
        for c in constlist:
            if c[0] not in beads_with_const:
                beads_with_const.append(c[0])
            if c[1] not in beads_with_const:
                beads_with_const.append(c[1])
        beads_with_const = sorted(beads_with_const)
        for i in range(len(beads_with_const)):
            for j in range(1 + i, len(beads_with_const)):
                const_exists = False
                for c in constlist:
                    if (c[0] == i and c[1] == j) or (c[0] == j and c[1] == i):
                        const_exists = True
                        break
                if not const_exists:
                    dist = np.linalg.norm(cgbead_coords[i] - cgbead_coords[j]) * 0.1
                    if dist < 0.35:
                        # Check that it's not in the bond list
                        in_bond_list = False
                        for b in bondlist:
                            if (b[0] == i and b[1] == j) or (b[0] == j and b[0] == i):
                                in_bond_list = True
                                break
                        # Are atoms part of the same ring
                        in_ring = False
                        for ring in ringatoms:
                            if cgbeads[i] in ring and cgbeads[j] in ring:
                                in_ring = True
                                break
                        # If not in bondlist and in the same ring, add the contraint
                        if not in_bond_list and in_ring:
                            constlist.append([i, j, dist])

        if not trial:

            if len(bondlist) > 0:
                text = "[bonds]\n" + "; i j     funct     length    force.c.\n"
                for b in bondlist:
                    # Make sure atoms in bond are not part of the same ring
                    text = text + "  {:d} {:d}       1         {:4.2f}     1250\n".format(
                        b[0] + 1, b[1] + 1, b[2]
                    )

            if len(constlist) > 0:
                text = text + "\n[constraints]\n" + ";  i   j     funct   length\n"

                for c in constlist:
                    text = text + "   {:<3d} {:<3d}   1       {:4.2f}\n".format(
                        c[0] + 1, c[1] + 1, c[2]
                    )
            # Make sure there's at least a bond to every atom
            for i in range(len(cgbeads)):
                bond_to_i = False
                for b in bondlist + constlist:
                    if i in [b[0], b[1]]:
                        bond_to_i = True
                if not bond_to_i:
                    print("Error. No bond to atom %d" % (i + 1))
                    exit(1)
    return bondlist, constlist, text


def print_angles(cgbeads, molecule, partitioning, cgbead_coords, bondlist, constlist, ringatoms):
    """print CG angles in itp format and returns the angles list"""
    logger.debug("Entering print_angles()")

    text = ""
    angle_list = []

    if len(cgbeads) > 2:
        # Angles
        for i in range(len(cgbeads)):
            for j in range(len(cgbeads)):
                for k in range(len(cgbeads)):
                    # Only up to 2 atoms can be ring-like
                    all_in_ring = False
                    for ring in ringatoms:
                        if cgbeads[i] in ring and cgbeads[j] in ring and cgbeads[k] in ring:
                            all_in_ring = True
                            break
                    # Forbid all atoms linked by constraints
                    all_constraints = False
                    ij_bonded = False
                    jk_bonded = False
                    ij_const = False
                    jk_const = False
                    for b in bondlist + constlist:
                        if i in [b[0], b[1]] and j in [b[0], b[1]]:
                            ij_bonded = True
                            if b in constlist:
                                ij_const = True
                        if j in [b[0], b[1]] and k in [b[0], b[1]]:
                            jk_bonded = True
                            if b in constlist:
                                jk_const = True
                    if ij_const and jk_const:
                        all_constraints = True
                    if (
                        not all_in_ring
                        and ij_bonded
                        and jk_bonded
                        and i != j
                        and j != k
                        and i != k
                        and not all_constraints
                    ):
                        # Measure angle between i, j, and k.
                        angle = (
                            180.0
                            / math.pi
                            * math.acos(
                                np.dot(
                                    cgbead_coords[i] - cgbead_coords[j],
                                    cgbead_coords[k] - cgbead_coords[j],
                                )
                                / (
                                    np.linalg.norm(cgbead_coords[i] - cgbead_coords[j])
                                    * np.linalg.norm(cgbead_coords[k] - cgbead_coords[j])
                                )
                            )
                        )
                        # Look for any double bond between atoms belonging to these CG beads.
                        atoms_in_fragment = []
                        for aa in partitioning.keys():
                            if partitioning[aa] == j:
                                atoms_in_fragment.append(aa)
                        forc_const = 25.0
                        for ib in range(len(molecule.GetBonds())):
                            abond = molecule.GetBondWithIdx(ib)
                            if (
                                abond.GetBeginAtomIdx() in atoms_in_fragment
                                and abond.GetEndAtomIdx() in atoms_in_fragment
                            ):
                                bondtype = molecule.GetBondBetweenAtoms(
                                    abond.GetBeginAtomIdx(), abond.GetEndAtomIdx()
                                ).GetBondType()
                                if bondtype == rdchem.BondType.DOUBLE:
                                    forc_const = 45.0
                        new_angle = True
                        for a in angle_list:
                            if i in a and j in a and k in a:
                                new_angle = False
                        if new_angle:
                            angle_list.append([i, j, k, angle, forc_const])

        if len(angle_list) > 0:
            text = text + "\n[angles]\n"
            text = text + "; i j k         funct   angle   force.c.\n"
            for a in angle_list:
                text = text + "  {:d} {:d} {:d}         2       {:<5.1f}  {:5.1f}\n".format(
                    a[0] + 1, a[1] + 1, a[2] + 1, a[3], a[4]
                )
            text = text + "\n"
    return text, angle_list


def print_dihedrals(cgbeads, constlist, ringatoms, cgbead_coords):
    """Print CG dihedrals in itp format"""
    logger.debug("Entering print_dihedrals()")

    text = ""

    if len(cgbeads) > 3:
        # Dihedrals
        dihed_list = []
        # Three ring atoms and one non ring
        for i in range(len(cgbeads)):
            for j in range(len(cgbeads)):
                for k in range(len(cgbeads)):
                    for l in range(len(cgbeads)):
                        if i != j and i != k and i != l and j != k and j != l and k != l:
                            # 3 atoms need to be ring like (in one ring!)
                            three_in_ring = False
                            for ring in ringatoms:
                                if [
                                    [cgbeads[i] in ring],
                                    [cgbeads[j] in ring],
                                    [cgbeads[k] in ring],
                                    [cgbeads[l] in ring],
                                ].count([True]) >= 3:
                                    three_in_ring = True
                                    break
                            for b in constlist:
                                if i in [b[0], b[1]] and j in [b[0], b[1]]:
                                    pass
                                if j in [b[0], b[1]] and k in [b[0], b[1]]:
                                    pass
                                if k in [b[0], b[1]] and l in [b[0], b[1]]:
                                    pass
                            # Distance criterion--beads can't be far apart
                            disthres = 0.35
                            close_enough = False
                            if (
                                np.linalg.norm(cgbead_coords[i] - cgbead_coords[j]) * 0.1
                                < disthres
                                and np.linalg.norm(cgbead_coords[j] - cgbead_coords[k]) * 0.1
                                < disthres
                                and np.linalg.norm(cgbead_coords[k] - cgbead_coords[l]) * 0.1
                                < disthres
                            ):
                                close_enough = True
                            already_dih = False
                            for dih in dihed_list:
                                if dih[0] == l and dih[1] == k and dih[2] == j and dih[3] == i:
                                    already_dih = True
                                    break
                            if three_in_ring and close_enough and not already_dih:
                                r1 = cgbead_coords[j] - cgbead_coords[i]
                                r2 = cgbead_coords[k] - cgbead_coords[j]
                                r3 = cgbead_coords[l] - cgbead_coords[k]
                                p1 = np.cross(r1, r2) / (np.linalg.norm(r1) * np.linalg.norm(r2))
                                p2 = np.cross(r2, r3) / (np.linalg.norm(r2) * np.linalg.norm(r3))
                                r2 /= np.linalg.norm(r2)
                                cosphi = np.dot(p1, p2)
                                sinphi = np.dot(r2, np.cross(p1, p2))
                                angle = 180.0 / math.pi * np.arctan2(sinphi, cosphi)
                                forc_const = 10.0
                                dihed_list.append([i, j, k, l, angle, forc_const])

        if len(dihed_list) > 0:
            text = text + "[dihedrals]"
            text = text + ";  i     j    k    l   funct   angle  force.c.\n"

            for d in dihed_list:
                text = (
                    text
                    + "   {:d}     {:d}    {:d}    {:d}       2     {:<5.1f}  {:5.1f}\n".format(
                        d[0] + 1, d[1] + 1, d[2] + 1, d[3] + 1, d[4], d[5]
                    )
                )
            text = text + "\n"

    return text


def smi2alogps(forcepred, smi, wc_log_p, bead, trial=False):
    """Returns water/octanol partitioning free energy
    according to ALOGPS"""
    logger.debug("Entering smi2alogps()")
    req = ""
    soup = ""
    try:
        session = requests.session()
        logger.debug("Calling http://vcclab.org/web/alogps/calc?SMILES=" + str(smi))
        req = session.get(
            "http://vcclab.org/web/alogps/calc?SMILES=" + str(smi.replace("#", "%23"))
        )
    except:
        print("Error. Can't reach vcclab.org to estimate free energy.")
        exit(1)
    try:
        doc = BeautifulSoup(req.content, "lxml")
    except Exception:
        raise
    try:
        soup = doc.prettify()
    except:
        print("Error with BeautifulSoup prettify")
        exit(1)
    found_mol_1 = False
    log_p = ""
    for line in soup.split("\n"):
        line = line.split()
        if "mol_1" in line:
            log_p = float(line[line.index("mol_1") + 1])
            found_mol_1 = True
            break
    if not found_mol_1:
        # If we're forcing a prediction, use Wildman-Crippen
        if forcepred:
            if trial:
                wrn = (
                    "; Warning: bead ID "
                    + str(bead)
                    + " predicted from Wildman-Crippen. Fragment "
                    + str(smi)
                    + "\n"
                )
                sys.stderr.write(wrn)
            log_p = wc_log_p
        else:
            print("ALOGPS can't predict fragment: %s" % smi)
            exit(1)
    logger.debug("logp value: %7.4f" % log_p)
    return convert_log_k(log_p)


def convert_log_k(log_k):
    """Convert log_{10}K to free energy (in kJ/mol)"""
    val = 0.008314 * 300.0 * log_k / math.log10(math.exp(1))
    logger.debug("free energy %7.4f kJ/mol" % val)
    return val


def mad(bead_type, delta_f, in_ring=False):
    """Mean absolute difference between type type and delta_f"""
    # logger.debug('Entering mad()')
    delta_f_types = read_delta_f_types()
    return math.fabs(delta_f_types[bead_type] - delta_f)


def determine_bead_type(delta_f, charge, hbonda, hbondd, in_ring):
    """Determine CG bead type from delta_f value, charge,
    and hbond acceptor, and donor"""
    if charge < -1 or charge > +1:
        print("Charge is too large: %s" % charge)
        print("No adequate force-field parameter.")
        exit(1)
    bead_type = []
    if in_ring:
        # We're in a ring, divide free energy by 3
        # (average number of beads per ring)
        if abs(delta_f) > 0.1:
            delta_f *= 2 / 3.0
    if charge != 0:
        # The compound has a +/- charge -> Q type
        if hbonda > 0 and hbondd > 0:
            bead_type = "Qda"
        elif hbonda > 0 and hbondd == 0:
            bead_type = "Qa"
        elif hbonda == 0 and hbondd > 0:
            bead_type = "Qd"
        else:
            bead_type = "Q0"
    else:
        # Neutral group
        # Use Hbond information only if we're close to Nda, Na, Nd types
        error = mad("Nda", delta_f, in_ring)
        if error < 3.0 and (hbonda > 0 or hbondd > 0):
            if hbonda > 0 and hbondd > 0:
                bead_type = "Nda"
            elif hbonda > 0 and hbondd == 0:
                bead_type = "Na"
            elif hbonda == 0 and hbondd > 0:
                bead_type = "Nd"
        else:
            # all other cases. Simply find the atom type that's closest in
            # free energy.
            othertypes = [
                "P5",
                "P4",
                "P3",
                "P2",
                "P1",
                "N0",
                "C5",
                "C4",
                "C3",
                "C2",
                "C1",
            ]
            min_error = 1000.0
            for cgtype in othertypes:
                tmp_error = mad(cgtype, delta_f, in_ring)
                if tmp_error < min_error:
                    min_error = tmp_error
                    bead_type = cgtype
            logger.debug("closest type: %s; error %7.4f" % (bead_type, min_error))
    if in_ring:
        bead_type = "S" + bead_type
    return bead_type

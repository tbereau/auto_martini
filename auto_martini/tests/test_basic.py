"""
Basic sanity test for the auto_martini package.
"""
import filecmp
import os
from pathlib import Path

import pytest

import auto_martini

dpath = Path("auto_martini/tests/files")


def test_auto_martini_imported():
    """Sample test, will always pass so long as import statement worked"""
    import sys

    assert "auto_martini" in sys.modules


@pytest.mark.parametrize(
    "sdf_file,num_beads", [(dpath / "benzene.sdf", 3), (dpath / "ibuprofen.sdf", 5)]
)
def test_auto_martini_run_sdf(sdf_file: str, num_beads: int):
    mol = auto_martini.topology.gen_molecule_sdf(str(sdf_file))
    cg_mol = auto_martini.solver.Cg_molecule(mol, "MOL")
    assert len(cg_mol.cg_bead_names) == num_beads


@pytest.mark.parametrize(
    "smiles,top_file,name,num_beads",
    [
        ("N1=C(N)NNC1N", "valid_GUA.top", "GUA", 2),
        ("CCC", "valid_PRO.top", "PRO", 1),
    ],
)
def test_auto_martini_run_smiles(smiles: str, top_file: Path, name: str, num_beads: int):
    mol, _ = auto_martini.topology.gen_molecule_smi(smiles)
    cg_mol = auto_martini.solver.Cg_molecule(mol, name)
    # assert filecmp.cmp(dpath / top_file, "mol.top")
    assert len(cg_mol.cg_bead_names) == num_beads

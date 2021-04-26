"""
Basic sanity test for the auto_martini package.
"""
import pytest
import auto_martini
import filecmp
from pathlib import Path
import os

dpath = Path("tests/files")

def test_auto_martini_imported():
    """Sample test, will always pass so long as import statement worked"""
    import sys

    assert "auto_martini" in sys.modules


@pytest.mark.parametrize("sdf_file", dpath.glob("*.sdf"))
def test_auto_martini_run_sdf(sdf_file):
    mol = auto_martini.topology.gen_molecule_sdf(str(sdf_file))
    auto_martini.solver.cg_molecule(mol, "MOL", "mol.top")
    Path("mol.top").unlink()


@pytest.mark.parametrize(
    "smiles,top_file,name",
    [
        ("N1=C(N)NNC1N", "valid_GUA.top", "GUA"),
        ("CCC", "valid_PRO.top", "PRO"),
    ]
)
def test_auto_martini_run_smiles(smiles: str, top_file: Path, name: str):
    mol, _ = auto_martini.topology.gen_molecule_smi(smiles)
    auto_martini.solver.cg_molecule(mol, name, "mol.top")
    # assert filecmp.cmp(dpath / top_file, "mol.top")
    # need to figure out a better way of comparing contents
    Path("mol.top").unlink()

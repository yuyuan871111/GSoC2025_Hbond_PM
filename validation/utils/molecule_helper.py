"""This script is copied and pasted from
https://github.com/chemosim-lab/ProLIF/pull/293.
"""

from collections.abc import Callable

from prolif.molecule import Molecule
from prolif.residue import ResidueId
from rdkit import Chem


def split_molecule(
    mol: Molecule, predicate: Callable[[ResidueId], bool]
) -> tuple[Molecule, Molecule]:
    """
    Splits a molecule into two based on a predicate function. The first molecule
    returned contains all residues of the input mol for which the predicate function was
    true, the second molecule contains the rest.
    This function is typically used to extract a ligand or water molecules from a file
    containing a solvated complex::
        >>> solvated_system = Molecule(...)
        >>> water_mol, protein_mol = split_molecule(
        ...     solvated_system, lambda x: x.name == "WAT"
        ... )
    .. versionadded:: 2.1.0
    """
    with Chem.RWMol(mol) as lhs, Chem.RWMol(mol) as rhs:
        for residue in mol:
            target = rhs if predicate(residue.resid) else lhs
            for atom in residue.GetAtoms():
                target.RemoveAtom(atom.GetUnsignedProp("mapindex"))
    return Molecule(lhs.GetMol()), Molecule(rhs.GetMol())

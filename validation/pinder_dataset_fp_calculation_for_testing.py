import pathlib

import dill
import pandas as pd
import prolif as plf
from prolif.io.protein_helper import ProteinHelper
from prolif.molecule import Molecule
from rdkit import Chem
from tqdm import tqdm

# validation set
with open("./data/pinder_test_local_paths.pkl", "rb") as f:
    local_paths = dill.load(f)

df = pd.read_csv("./data/pinder_test_systems.csv")

protein_helper = ProteinHelper(
    [
        {
            "MSE": {"SMILES": "C[Se]CC[CH](N)C=O"},
            "HIC": {"SMILES": "Cn1cc(nc1)C[C@@H](C(=O))N"},
            "OAS": {"SMILES": "CC(=O)OC[C@@H](C(=O))N"},
            "KPI": {"SMILES": r"C/C(=N\\CCCC[C@@H](C(=O))N)/C(=O)O"},
            "MIS": {"SMILES": "CC(C)O[P@](=O)(O)OC[C@@H](C(=O))N"},
            "TPO": {"SMILES": "C[C@H]([C@@H](C(=O))N)OP(=O)(O)O"},
            "SEP": {"SMILES": "C([C@@H](C(=O))N)OP(=O)(O)O"},
            "PYR": {"SMILES": "CC(=O)C(=O)"},
            "SNN": {"SMILES": "O=CC(N)CC(=O)"},  # prune the N in the ring
            "FME": {"SMILES": "CSCC[C@@H](C(=O))NC=O"},
            "KCX": {"SMILES": "C(CCNC(=O)O)C[C@@H](C(=O))N"},
            "MHS": {"SMILES": "Cn1cncc1C[C@@H](C(=O))N"},
            "SME": {"SMILES": "C[S@@](=O)CC[C@@H](C(=O))N"},
            "PCA": {"SMILES": "C1CC(=O)N[C@@H]1C(=O)"},
            "ACE": {"SMILES": "CC=O"},
            "MHI": {"SMILES": "c1cc(ccc1C[C@@H](C(=O))N)I"},
            "9IJ": {"SMILES": "c1ccc(c(c1)C[C@@H](C(=O))N)C#N"},
            "CGU": {"SMILES": "C(C(C(=O)O)C(=O)O)[C@@H](C(=O))N"},
            "CSO": {"SMILES": "C([C@@H](C(=O))N)SO"},
            "ALY": {"SMILES": "CC(=O)NCCCC[C@@H](C(=O))N"},
            "PHD": {"SMILES": "C([C@@H](C(=O))N)C(=O)OP(=O)(O)O"},
            "OCS": {"SMILES": "C([C@@H](C(=O))N)S(=O)(=O)O"},
            "MDO": {"SMILES": "C[C@@H](C1=NC(=C)C(=O)N1CC(=O))N"},
            "SAR": {"SMILES": "CNCC(=O)"},
            "AGM": {"SMILES": "C[C@@H](CC[C@@H](C(=O))N)NC(=[NH2+])N"},
            "PHI": {"SMILES": "c1cc(ccc1C[C@@H](C(=O))N)I"},
            "0AF": {"SMILES": "c1cc2c(c[nH]c2c(c1)O)C[C@@H](C(=O))N"},
            "MVA": {"SMILES": "CC(C)[C@@H](C(=O))NC"},
            "SMC": {"SMILES": "CSC[C@@H](C(=O))N"},
            "GL3": {"SMILES": "C(C(=O)S)N"},
            "IML": {"SMILES": "CC[C@H](C)[C@@H](C(=O))NC"},
            "I2M": {"SMILES": "CCC(C)(C)[C@@H](C(=O))N"},
            "MGN": {"SMILES": "C[C@](CCC(=O)N)(C(=O))N"},
        },
    ]
)


for idx in tqdm(df.index):
    each_system_id = df.loc[idx].id
    each_system = local_paths[each_system_id]

    raw_receptor_path = each_system["holo_receptor"]
    receptor_path = str(
        raw_receptor_path.parent.parent
        / f"protonated_{raw_receptor_path.parent.stem}"
        / f"{raw_receptor_path.stem}_protonated.pdb"
    )
    raw_ligand_path = each_system["holo_ligand"]
    ligand_path = str(
        raw_ligand_path.parent.parent
        / f"protonated_{raw_ligand_path.parent.stem}"
        / f"{raw_ligand_path.stem}_protonated.pdb"
    )

    # Check if the fingerprint files already exist
    if (
        pathlib.Path(f"./pinder_test/explicit/fp_{idx}.pkl").exists()
        and pathlib.Path(f"./pinder_test/implicit/fp_{idx}.pkl").exists()
    ):
        continue

    # Load files
    try:
        protein_mol = protein_helper.standardize_protein(
            Molecule.from_rdkit(
                Chem.MolFromPDBFile(
                    receptor_path,
                    sanitize=False,
                    removeHs=False,
                    proximityBonding=True,
                )
            )
        )
        ligand = protein_helper.standardize_protein(
            Molecule.from_rdkit(
                Chem.MolFromPDBFile(
                    ligand_path,
                    sanitize=False,
                    removeHs=False,
                    proximityBonding=True,
                )
            )
        )

    except Exception as e:
        print(f"*** Error processing {each_system_id}: {e}")  # noqa: T201
        continue

    try:
        # [explicit hbond calculation]
        # Calculate explicit hydrogen bonds
        if not pathlib.Path(f"./pinder_test/explicit/fp_{idx}.pkl").exists():
            fp = plf.Fingerprint(["HBDonor", "HBAcceptor"], count=True)
            fp.run_from_iterable([ligand], protein_mol, progress=False)
            fp.to_pickle(f"./pinder_test/explicit/fp_{idx}.pkl")

        # [implicit hbond calculation]
        # Remove hydrogens from the protein and ligand
        ligand_i = protein_helper.standardize_protein(
            Molecule.from_rdkit(Chem.RemoveAllHs(ligand, sanitize=False))
        )
        protein_mol_i = protein_helper.standardize_protein(
            Molecule.from_rdkit(Chem.RemoveAllHs(protein_mol, sanitize=False))
        )

        # Calculate implicit hydrogen bonds
        if not pathlib.Path(f"./pinder_test/implicit/fp_{idx}.pkl").exists():
            fp = plf.Fingerprint(
                ["ImplicitHBDonor", "ImplicitHBAcceptor"],
                count=True,
            )
            fp.run_from_iterable([ligand_i], protein_mol_i, progress=False)
            fp.to_pickle(f"./pinder_test/implicit/fp_{idx}.pkl")

    except Exception as e:
        print("DEBUG later:")  # noqa: T201
        print(f"*** Error processing {each_system_id}: {e}")  # noqa: T201
        continue

import pathlib

import dill
import pandas as pd
import prolif as plf
from prolif.io.protein_helper import ProteinHelper
from prolif.molecule import Molecule
from rdkit import Chem
from tqdm import tqdm
from utils.molecule_helper import split_molecule

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
            "CME": {"SMILES": "C(CSSC[C@@H](C(=O))N)O"},
            "MLZ": {"SMILES": "CNCCCC[C@@H](C(=O))N"},
            "LLP": {"SMILES": "Cc1c(c(c(cn1)COP(=O)(O)O)/C=N/CCCCC(C(=O))N)O"},
            "DAL": {"SMILES": "C[C@H](C(=O))N"},
            "DDZ": {"SMILES": "[C@H](C(O)O)(C(=O))N"},
            "LP6": {"SMILES": "C1CCN(CC1)CCCC[C@@H](C(=O))N"},
            "CSX": {"SMILES": "C([C@@H](C(=O))N)[S@H]=O"},
            "NEP": {"SMILES": "c1c(ncn1P(=O)(O)O)C[C@@H](C(=O))N"},
            "5VV": {"SMILES": "C([C@@H](C(=O)O)NC(=O))C(=O)N"},
            "2CO": {"SMILES": "C([C@@H](C(=O))N)SOO"},
            "CSA": {"SMILES": "CC(=O)CSC[C@@H](C(=O))N"},
            "MLY": {"SMILES": "CN(C)CCCC[C@@H](C(=O))N"},
            "LYR": {
                "SMILES": (
                    "CC1=C(C(CCC1)(C)C)"
                    r"C=C/C(=C/C=C\\C(=C\\CNCCCC[C@@H](C(=O))N)\\C)/C"
                )
            },
            "FGP": {"SMILES": "[C@H]([C@@H](O)OP(=O)(O)O)(C(=O))N"},
            "CAS": {"SMILES": "C[As](C)SC[C@@H](C(=O))N"},
            "SCH": {"SMILES": "CSSC[C@@H](C(=O))N"},
            "V44": {"SMILES": "Cn1ccnc1C(CSC[C@@H](C(=O))N)c2nccn2C"},
            "SEB": {"SMILES": "c1ccc(cc1)CS(=O)(=O)OC[C@@H](C(=O))N"},
            "IYR": {"SMILES": "c1cc(c(cc1C[C@@H](C(=O))N)I)O"},
            "OHI": {"SMILES": "C1=NC(=O)N=C1C[C@@H](C(=O))N"},
            "CSS": {"SMILES": "C([C@@H](C(=O))N)SS"},
            "CR2": {
                "SMILES": r"c1cc(ccc1\\C=C/2\\C(=O)N(C(=N2)CN)CC(=O))O"
            },  # need check
            "TYI": {"SMILES": "c1c(cc(c(c1I)O)I)C[C@@H](C(=O))N"},
            "CSD": {"SMILES": "C([C@@H](C(=O))N)[S@@](=O)O"},
            "MSO": {"SMILES": "C[Se@@](=O)CC[C@@H](C(=O))N"},
            "2ML": {"SMILES": "CC(C)C[C@@](C)(C(=O))N"},
            "R1A": {"SMILES": "CC1(C=C(C([N+]1=O)(C)C)CSSC[C@@H](C(=O))N)C"},
            "PTR": {"SMILES": "c1cc(ccc1C[C@@H](C(=O))N)OP(=O)(O)O"},
            "4MM": {"SMILES": "C[N+](C)(C)[C@@H](CCSC)C(=O)"},
            "IO8": {"SMILES": r"c1ccc2c(c1)c(c[nH]2)/C=C\\3/C(=O)N(C(=N3)CN)CC(=O)"},
            "M0H": {"SMILES": "C([C@@H](C(=O))N)SCO"},
            "SUN": {"SMILES": "CCO[P@](=O)(N(C)C)OC[C@@H](C(=O))N"},
            "AYA": {"SMILES": "C[C@@H](C(=O))NC(=O)C"},
            "0YG": {"SMILES": r"c1cc(ccc1/C=C(/C(=O)NCC(=O))\\N)O"},
            "FTR": {"SMILES": "c1cc2c(cc1F)c(c[nH]2)C[C@@H](C(=O))N"},
            "YCM": {"SMILES": "C([C@@H](C(=O))N)SCC(=O)N"},
            "AGD": {"SMILES": "c1nc2c(n1C[C@H](C(=O))N)N=C(NC2=O)N"},
            "DAH": {"SMILES": "c1cc(c(cc1C[C@@H](C(=O))N)O)O"},
            "MEQ": {"SMILES": "CNC(=O)CC[C@@H](C(=O))N"},
            "HTR": {"SMILES": "c1ccc2c(c1)c(c[nH]2)[C@@H]([C@@H](C(=O))N)O"},
            "PBF": {"SMILES": "c1ccc(cc1)C(=O)c2ccc(cc2)C[C@@H](C(=O))N"},
            "SCY": {"SMILES": "CC(=O)SC[C@@H](C(=O))N"},
            "TYQ": {"SMILES": "	c1c(c(cc(c1N)O)O)C[C@@H](C(=O))N"},
        },
    ]
)


for idx in tqdm(df.index):
    each_system_id = df.loc[idx].id
    each_system = local_paths[each_system_id]

    # raw_receptor_path = each_system["holo_receptor"]
    # receptor_path = str(
    #     raw_receptor_path.parent.parent
    #     / f"protonated_{raw_receptor_path.parent.stem}"
    #     / f"{raw_receptor_path.stem}_protonated.pdb"
    # )
    # raw_ligand_path = each_system["holo_ligand"]
    # ligand_path = str(
    #     raw_ligand_path.parent.parent
    #     / f"protonated_{raw_ligand_path.parent.stem}"
    #     / f"{raw_ligand_path.stem}_protonated.pdb"
    # )
    pdbs_dir = each_system["holo_receptor"].parent.parent / "protonated_pdbs"
    native_pdb_path = pdbs_dir / f"{each_system_id}_protonated.pdb"

    # Check if the fingerprint files already exist
    if (
        pathlib.Path(f"./pinder_test/explicit/fp_{idx}.pkl").exists()
        and pathlib.Path(f"./pinder_test/implicit/fp_{idx}.pkl").exists()
    ):
        continue

    # Load files
    try:
        pl_complex = Molecule.from_rdkit(
            Chem.MolFromPDBFile(
                native_pdb_path,
                sanitize=False,
                removeHs=False,
                proximityBonding=True,
            )
        )
        ligand, receptor = split_molecule(pl_complex, lambda x: x.chain == "L")
        ligand, _ = split_molecule(ligand, lambda x: x.chain == "L")
        receptor, _ = split_molecule(receptor, lambda x: x.chain == "R")
        receptor = protein_helper.standardize_protein(receptor)
        ligand = protein_helper.standardize_protein(ligand)

    except Exception as e:
        print(f"*** Error processing {each_system_id}: {e}")  # noqa: T201
        continue

    try:
        # [explicit hbond calculation]
        # Calculate explicit hydrogen bonds
        if not pathlib.Path(f"./pinder_test/explicit/fp_{idx}.pkl").exists():
            fp = plf.Fingerprint(["HBDonor", "HBAcceptor"], count=True)
            fp.run_from_iterable([ligand], receptor, progress=False)
            fp.to_pickle(f"./pinder_test/explicit/fp_{idx}.pkl")

        # [implicit hbond calculation]
        # Remove hydrogens from the protein and ligand
        ligand_i = protein_helper.standardize_protein(
            Molecule.from_rdkit(Chem.RemoveAllHs(ligand, sanitize=False))
        )
        receptor_i = protein_helper.standardize_protein(
            Molecule.from_rdkit(Chem.RemoveAllHs(receptor, sanitize=False))
        )

        # Calculate implicit hydrogen bonds
        if not pathlib.Path(f"./pinder_test/implicit/fp_{idx}.pkl").exists():
            fp = plf.Fingerprint(
                ["ImplicitHBDonor", "ImplicitHBAcceptor"],
                count=True,
            )
            fp.run_from_iterable([ligand_i], receptor_i, progress=False)
            fp.to_pickle(f"./pinder_test/implicit/fp_{idx}.pkl")

    except Exception as e:
        print("DEBUG later:")  # noqa: T201
        print(f"*** Error processing {each_system_id}: {e}")  # noqa: T201
        continue

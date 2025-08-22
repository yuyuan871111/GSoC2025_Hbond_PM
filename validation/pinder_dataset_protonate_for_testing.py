import dill
from tqdm import tqdm
from utils.protonate import pinder_system_protonate

# validation set
with open("./data/pinder_test_local_paths.pkl", "rb") as f:
    local_paths = dill.load(f)

for each_system_id in tqdm(local_paths):
    each_system = local_paths[each_system_id]

    # receptor
    # receptor_path = each_system["holo_receptor"]
    # pinder_system_protonate(input_pdb=receptor_path)

    # ligand
    # ligand_path = each_system["holo_ligand"]
    # pinder_system_protonate(input_pdb=ligand_path)

    # receptor and ligand in test set is aligned (cannot find the non-aligned one)
    # instead, use the native pdb files
    pdbs_dir = each_system["holo_receptor"].parent.parent / "pdbs"
    native_pdb_path = pdbs_dir / f"{each_system_id}.pdb"
    pinder_system_protonate(input_pdb=native_pdb_path)

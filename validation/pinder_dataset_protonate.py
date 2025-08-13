import dill
from tqdm import tqdm
from utils.protonate import pinder_system_protonate

# validation set
with open("./data/pinder_val_local_paths.pkl", "rb") as f:
    local_paths = dill.load(f)

for each_system_id in tqdm(local_paths):
    each_system = local_paths[each_system_id]

    # receptor
    receptor_path = each_system["holo_receptor"]
    pinder_system_protonate(input_pdb=receptor_path)

    # ligand
    ligand_path = each_system["holo_ligand"]
    pinder_system_protonate(input_pdb=ligand_path)

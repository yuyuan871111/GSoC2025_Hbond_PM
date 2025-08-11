import pandas as pd
from tqdm import tqdm
from utils.protonate import plinder_system_protonate

df = pd.read_csv("./data/plinder_test_systems.csv")

# each_system_id = df.iloc[0]['system_id']  # Example to access system_id
for each_system_id in tqdm(df["system_id"].unique()):
    # Create a PlinderSystem instance for each system_id
    # plinder_system = PlinderSystem(system_id=each_system_id)
    plinder_system_protonate(plinder_system_id=each_system_id)

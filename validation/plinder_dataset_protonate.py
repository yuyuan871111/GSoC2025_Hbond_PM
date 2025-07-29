
import pandas as pd

from tqdm import tqdm
from protonate import plinder_system_protonate

df = pd.read_csv("plinder_val_systems.csv")

# each_system_id = df.iloc[0]['system_id']  # Example to access system_id
for idx, each_system_id in enumerate(tqdm(df['system_id'])):

    # Create a PlinderSystem instance for each system_id
    # plinder_system = PlinderSystem(system_id=each_system_id)
    plinder_system_protonate(plinder_system_id=each_system_id)
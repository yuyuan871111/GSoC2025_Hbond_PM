# global_vars.py: Global variables used by MaSIF --
# mainly pointing to environment variables of programs used by MaSIF.
# Pablo Gainza - LPDI STI EPFL 2018-2019
# Released under an Apache License 2.0
# Modified by Yu-Yuan Yang (2025)

import configparser

# Read config file.
config = configparser.ConfigParser()
path_args = __file__.split("/")[0:-1]
root_path = "/".join(path_args)
config.read(f"{root_path}/config.cfg")  # serve as a main module
config.sections()

# Set the environment variables for the programs used by MaSIF.
reduce_bin = ""
if "REDUCE_BIN" in config["ThirdParty"]:
    reduce_bin = config["ThirdParty"]["REDUCE_BIN"]
else:
    print("Warning: REDUCE_BIN not set. Variable should point to REDUCE program.")  # noqa: T201

pdb2pqr_bin = ""
if "PDB2PQR_BIN" in config["ThirdParty"]:
    pdb2pqr_bin = config["ThirdParty"]["PDB2PQR_BIN"]
else:
    print("Warning: PDB2PQR_BIN not set. Variable should point to PDB2PQR_BIN program.")  # noqa: T201

class NoSolutionError(Exception):
    pass

"""
protonate.py: Wrapper method for the reduce program: protonate (i.e., add hydrogens) 
a pdb using reduce (or propka) and save to an output file.
Pablo Gainza - LPDI STI EPFL 2019
Released under an Apache License 2.0

Modified by Yu-Yuan Yang (2025)
"""
import os
import random
from subprocess import PIPE, Popen
from time import strftime

import plinder.core.utils.config
from global_vars import pdb2pqr_bin, reduce_bin

cfg = plinder.core.get_config()


def protonate(
        in_pdb_file: str, 
        out_pdb_file: str, 
        method: str = "propka", 
        keep_tempfiles: bool = False,
        ):
    """
    Protonate (i.e., add hydrogens) a pdb using reduce or propka and
    save to an output file.
    
    Args
    ----
    in_pdb_file: file to protonate.
    out_pdb_file: output file where to save the protonated pdb file.
    method: method to use for protonation. Options are "reduce" or "propka".
            If None, no protonation is performed.
    """
    if method is None:
        pass
    elif method == "reduce":
        # Remove protons first, in case the structure is already protonated
        args = [reduce_bin, "-Trim", in_pdb_file]
        p2 = Popen(args, stdout=PIPE, stderr=PIPE)
        stdout, stderr = p2.communicate()
        outfile = open(out_pdb_file, "w")
        outfile.write(stdout.decode("utf-8").rstrip())
        outfile.close()
        # Now add them again.
        args = [reduce_bin, "-HIS", out_pdb_file]
        p2 = Popen(args, stdout=PIPE, stderr=PIPE)
        stdout, stderr = p2.communicate()
        # write the output to a file.
        outfile = open(out_pdb_file, "w")
        outfile.write(stdout.decode("utf-8"))
        outfile.close()

    elif method == "propka":
        # Create temporary directory
        now = strftime("%y%m%d%H%M%S")
        randnum = str(random.randint(1, 10000000))
        tmp_file_base = f"{in_pdb_file.replace('.pdb', '')}_temp_{now}_{randnum}"
        os.makedirs(tmp_file_base, exist_ok=False)

        filename = out_pdb_file.replace(".pdb", "").split("/")[-1]
        apbs_in_file = f"{tmp_file_base}/{filename}.in"
        pqr_out_file = f"{tmp_file_base}/{filename}"

        args = [
            pdb2pqr_bin,
            "--ff=PARSE",
            "--ffout=AMBER",
            "--pdb-output",
            out_pdb_file,
            "--whitespace",
            "--apbs-input",
            apbs_in_file,
            "--titration-state-method=propka",
            in_pdb_file,
            pqr_out_file,
        ]
        p2 = Popen(args, stdout=PIPE, stderr=PIPE)
        stdout, stderr = p2.communicate()

        if os.path.exists(out_pdb_file):
            if keep_tempfiles:
                return tmp_file_base
            else:
                os.system(f"rm -r {tmp_file_base}")
        else:
            raise Exception(f"Error in protonation. Check log file in {tmp_file_base}")

    else:
        raise ValueError(f"Unknown protonation method: {method}")


def protonate_ligand(in_sdf_file, out_sdf_file):
    """
    Protonate a ligand in sdf format using and save to an output file.
    in_sdf_file: file to protonate.
    out_sdf_file: output file where to save the protonated sdf file.
    """
    args = ["obabel", "-isdf", in_sdf_file, "-osdf", "-O", out_sdf_file, "-h"]
    p2 = Popen(args, stdout=PIPE, stderr=PIPE)
    stdout, stderr = p2.communicate()

    if p2.returncode != 0:
        
        with open(f"{out_sdf_file}.log", "w") as log_file:
            log_file.write(stdout.decode("utf-8"))
            log_file.write(stderr.decode("utf-8"))
        
        raise Exception(f"Error in protonation. Check log file in {out_sdf_file}.log")


def plinder_system_protonate(plinder_system_id):
    # I/O
    # input_pdb = plinder_system.receptor_pdb
    input_pdb = f"{cfg.data.plinder_dir}/systems/{plinder_system_id}/receptor.pdb"

    # output_folder = os.path.dirname(input_pdb).replace("systems", "protonated_systems")
    output_folder = f"{cfg.data.plinder_dir}/protonated_systems/{plinder_system_id}"
    os.makedirs(output_folder, exist_ok=True)
    output_pdb = f"{output_folder}/receptor_protonated.pdb"

    # protonate for the receptor
    protonate(
        in_pdb_file=input_pdb,
        out_pdb_file=output_pdb,
        method="propka"
    )

    # protonate for the ligands
    ligand_file_path = f"{cfg.data.plinder_dir}/systems/{plinder_system_id}/ligand_files"
    ligand_sdfs = os.listdir(ligand_file_path)
    for ligand_sdf in ligand_sdfs:

        ligand_sdf_path = f"{ligand_file_path}/{ligand_sdf}"
        output_ligand_sdf = f"{output_folder}/{ligand_sdf.replace('.sdf', '_protonated.sdf')}"
        
        protonate_ligand(
            in_sdf_file=ligand_sdf_path,
            out_sdf_file=output_ligand_sdf
        )

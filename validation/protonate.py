"""
protonate.py: Wrapper method for the reduce program: protonate (i.e., add hydrogens)
a pdb using reduce (or propka) and save to an output file.
Pablo Gainza - LPDI STI EPFL 2019
Released under an Apache License 2.0

Modified by Yu-Yuan Yang (2025)
"""

import os
import pathlib
import random
from subprocess import PIPE, Popen
from time import strftime

import MDAnalysis as mda
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
        with open(out_pdb_file, "w") as f:
            f.write(stdout.decode("utf-8").rstrip())
        # Now add them again.
        args = [reduce_bin, "-HIS", out_pdb_file]
        p2 = Popen(args, stdout=PIPE, stderr=PIPE)
        stdout, stderr = p2.communicate()
        # write the output to a file.
        with open(out_pdb_file, "w") as f:
            f.write(stdout.decode("utf-8"))

    elif method == "propka":
        # Create temporary directory
        now = strftime("%y%m%d%H%M%S")
        randnum = str(random.randint(1, 10000000))
        tmp_file_base = f"{in_pdb_file.replace('.pdb', '')}_temp_{now}_{randnum}"
        os.makedirs(tmp_file_base, exist_ok=False)  # noqa: PTH103

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
        stdout, stderr = p2.communicate()  # noqa: RUF059, RUF100

        if pathlib.Path(out_pdb_file).exists() and p2.returncode == 0:
            if keep_tempfiles:
                return
            os.system(f"rm -r {tmp_file_base}")
        else:
            raise Exception(f"Error in protonation. Check log file in {tmp_file_base}")

    else:
        raise ValueError(f"Unknown protonation method: {method}")


def protonate_with_propka_without_opt(
    in_pdb_file: str,
    out_pdb_file: str,
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

    # Create temporary directory
    now = strftime("%y%m%d%H%M%S")
    randnum = str(random.randint(1, 10000000))
    tmp_file_base = f"{in_pdb_file.replace('.pdb', '')}_temp_{now}_{randnum}"
    os.makedirs(tmp_file_base, exist_ok=False)  # noqa: PTH103

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
        "--noopt",
        "--apbs-input",
        apbs_in_file,
        "--titration-state-method=propka",
        in_pdb_file,
        pqr_out_file,
    ]
    p2 = Popen(args, stdout=PIPE, stderr=PIPE)
    stdout, stderr = p2.communicate()  # noqa: RUF059, RUF100

    if pathlib.Path(out_pdb_file).exists():
        if keep_tempfiles:
            return
        os.system(f"rm -r {tmp_file_base}")
    else:
        raise Exception(f"Error in protonation. Check log file in {tmp_file_base}")


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


def check_errors(log_file):
    with open(log_file) as file:
        lines = file.readlines()

        return {
            "log_file": log_file,
            "raw_error_msg": lines,
            "error_msg": [line.strip() for line in lines if "CRITICAL:" in line],
        }


def try_protonate_with_reduce(
    input_pdb: str,
    output_pdb: str,
):
    """
    Try to protonate a system with reduce.
    If it fails, it will try to protonate with propka.
    """
    try:
        print("Trying to protonate with reduce.")  # noqa: T201
        protonate(in_pdb_file=input_pdb, out_pdb_file=output_pdb, method="reduce")
        print("Protonation with reduce was successful.")  # noqa: T201
        return

    except Exception as e:
        print(f"[TODO] Protonation with reduce failed: {e}")  # noqa: T201
        return


def error_handler(errors: dict, output_pdb: str | None = None):
    plinder_system_id = errors["log_file"].split("/")[-3]
    input_pdb = f"{cfg.data.plinder_dir}/systems/{plinder_system_id}/receptor.pdb"

    # If not the error from PDB2PQR (or said from PROPKA)
    # Try to protonate with reduce
    if errors.get("error_msg", []) == []:
        try_protonate_with_reduce(
            input_pdb=input_pdb,
            output_pdb=output_pdb,
        )

    # If the error is from PDB2PQR
    # Check the error message and handle accordingly
    for error_msg in errors.get("error_msg", []):
        if "Too few atoms present to reconstruct or cap residue" in error_msg:
            print(  # noqa: T201
                "Trimming the receptor to remove residue with too many missing "
                "heavy atoms."
            )
            temp = error_msg.split("Heavy atoms missing from ")[1]
            residue = temp.split(":")[0].split(" ")
            resname = residue[0]
            chain = residue[1]
            resid = residue[2]

            u = mda.Universe(input_pdb)

            trimed_file_path = (
                f"{cfg.data.plinder_dir}/systems/{plinder_system_id}/receptor_trim.pdb"
            )
            u.select_atoms(
                f"not (resname {resname} and segid {chain} and resid {resid})"
            ).write(trimed_file_path)

            print(  # noqa: T201
                f"Trimmed {resname} {chain} {resid} "
                f"from {plinder_system_id}/receptor.pdb."
            )

            try:
                print("Trying to protonate with propka after trimming.")  # noqa: T201
                protonate(
                    in_pdb_file=trimed_file_path,
                    out_pdb_file=output_pdb,
                    method="propka",
                )
                print("Protonation with propka was successful after trimming.")  # noqa: T201
                return

            except Exception as e:
                print(f"Protonation with propka after trimming failed: {e}")  # noqa: T201
                try_protonate_with_reduce(
                    input_pdb=trimed_file_path,
                    output_pdb=output_pdb,
                )
                return

        elif (
            "cannot convert float NaN to integer" in error_msg
            or "from integral, exceeding error tolerance" in error_msg
            or "Biomolecular structure is incomplete" in error_msg
        ):
            try:
                print("Trying to protonate with propka without optimization.")  # noqa: T201
                protonate_with_propka_without_opt(
                    in_pdb_file=input_pdb, out_pdb_file=output_pdb
                )
                print("Protonation with propka without optimization was successful.")  # noqa: T201
                return

            except Exception as e:
                print(f"Protonation with propka without optimization failed: {e}")  # noqa: T201
                try_protonate_with_reduce(
                    input_pdb=input_pdb,
                    output_pdb=output_pdb,
                )
                return

        else:
            print(f"[TODO] Cannot handle this error yet: {error_msg}.")  # noqa: T201
            return


def plinder_system_protonate(plinder_system_id, p_method="propka"):
    # I/O
    # input_pdb = plinder_system.receptor_pdb
    input_pdb = f"{cfg.data.plinder_dir}/systems/{plinder_system_id}/receptor.pdb"

    # output_folder = os.path.dirname(input_pdb).replace(
    # "systems", "protonated_systems")
    output_folder = f"{cfg.data.plinder_dir}/protonated_systems/{plinder_system_id}"
    os.makedirs(output_folder, exist_ok=True)  # noqa: PTH103
    output_pdb = f"{output_folder}/receptor_protonated.pdb"

    # protonate for the receptor
    if not pathlib.Path(output_pdb).exists():
        try:
            protonate(in_pdb_file=input_pdb, out_pdb_file=output_pdb, method=p_method)
        except Exception as e:
            print(f"Error occurred while protonating receptor: {e}")  # noqa: T201
            e_log_file = (
                f"{e.args[0].split('Check log file in ')[-1].strip()}/"
                "receptor_protonated.log"
            )
            errors = check_errors(e_log_file)
            error_handler(errors, output_pdb=output_pdb)

    # protonate for the ligands
    ligand_file_path = (
        f"{cfg.data.plinder_dir}/systems/{plinder_system_id}/ligand_files"
    )
    ligand_sdfs = os.listdir(ligand_file_path)  # noqa: PTH208
    for ligand_sdf in ligand_sdfs:
        ligand_sdf_path = f"{ligand_file_path}/{ligand_sdf}"
        output_ligand_sdf = (
            f"{output_folder}/{ligand_sdf.replace('.sdf', '_protonated.sdf')}"
        )

        if not pathlib.Path(output_ligand_sdf).exists():
            protonate_ligand(
                in_sdf_file=ligand_sdf_path, out_sdf_file=output_ligand_sdf
            )

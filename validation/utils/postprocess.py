import re
from typing import cast

import numpy as np
import prolif as plf
from prolif.residue import ResidueId

_RE_RESID = re.compile(
    r"(TIP[234]|T[234]P|H2O|[0-9][A-Z]{2}|[A-Z ]+)?(\d*)\.?([A-Z_\d]{1,2})?"
)


def residue_id_from_string(resid_str: str) -> ResidueId:
    matches = cast(re.Match, _RE_RESID.search(resid_str))
    name, number, chain = matches.groups()
    number = int(number) if number else 0
    return ResidueId(name, number, chain)


def get_metadata_from_implicit_hbond_using_set(
    a_set: set,
    fp_i: plf.Fingerprint,
    ignore_water: bool = True,
):
    aaa_dev_list = []
    daa_dev_list = []
    dpa_list = []
    apa_list = []
    vina_hbond_potential_list = []
    for each_pair in a_set:
        first_residue = residue_id_from_string(each_pair[0])
        second_residue = residue_id_from_string(each_pair[1])
        if ignore_water and "HOH" in {first_residue.name, second_residue.name}:
            continue
        metadata = fp_i.ifp[0].get((first_residue, second_residue), {})
        for each_interaction_type in metadata:
            best_aaa_dev = np.inf
            best_daa_dev = np.inf
            best_dev_square_sum = np.inf
            best_vina_hbond_potential = np.inf
            for each_metadata in metadata[each_interaction_type]:
                aaa_dev = each_metadata["acceptor_atom_angle_deviation"]
                daa_dev = each_metadata["donor_atom_angle_deviation"]
                vina_hbond_potential = each_metadata["vina_hbond_potential"]
                dev_square_sum = aaa_dev**2 + daa_dev**2

                if dev_square_sum < best_dev_square_sum:
                    best_dev_square_sum = dev_square_sum
                    best_aaa_dev = aaa_dev
                    best_daa_dev = daa_dev
                    best_vina_hbond_potential = vina_hbond_potential

                if "donor_plane_angle" in each_metadata:
                    dpa_list.append(each_metadata["donor_plane_angle"])
                if "acceptor_plane_angle" in each_metadata:
                    apa_list.append(each_metadata["acceptor_plane_angle"])

            # only the least deviated angles are recorded
            aaa_dev_list.append(best_aaa_dev)
            daa_dev_list.append(best_daa_dev)
            vina_hbond_potential_list.append(best_vina_hbond_potential)

    return (aaa_dev_list, daa_dev_list, dpa_list, apa_list, vina_hbond_potential_list)

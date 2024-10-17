import numpy as np
from scipy.stats import rankdata
from molbar.barcodes.barcodes import _calculate_barcode
from molbar.molgraph.nodes.priorities import Priorities

from molbar.indices.absolute_conf import AbsoluteConfIndex


def determine_absolute_configuration(molgraph):
    """
    Iterate again over all fragments and calculate the absolute configuration index to account for pseudochirality.
    Args:
        fragments_data (dict): container for all fragment data

    Returns:
        dict: container for all fragment data with updated absolute configuration index
    """
    fragments_data = molgraph.fragments_data
    molgraph.set_nodes_visible(include_all=True)
    atomic_numbers = molgraph.return_node_data(attribute="atomic_numbers").astype(float)
    for fragment_id, fragment_data in fragments_data.items():
        molgraph.set_nodes_visible(visible_nodes=fragment_data["visible_nodes"])
        geometry = fragment_data["final_geometry"]
        priorities = molgraph.return_node_data(attribute="priorities")
        core_indices = fragment_data["core_indices"]
        absconf = AbsoluteConfIndex(geometry, priorities, core_indices)
        G0s = absconf.get_absconf_index()
        fragment_data["absconf_index"] = G0s
        if G0s != 0:
            if G0s == 1:
                shift = 0.49
            else:
                shift = -0.49
            global_core_indices = np.array(fragment_data["visible_nodes"])[core_indices]
            atomic_numbers[global_core_indices] += shift
    molgraph.set_nodes_visible(include_all=True)
    n_atoms = molgraph.return_n_atoms()
    cm_matrix = np.round(molgraph.unified_coulomb_matrix, 4)
    degrees = molgraph.get_degree()
    eps, C = _calculate_barcode(
        n_atoms, cm_matrix, atomic_numbers, degrees, return_eigenvec=True
    )
    densities = Priorities().get_densities(eps, C, 5)
    molgraph.set_nodes_visible(include_all=True)
    priorities = rankdata(densities, method="min")
    for fragment_id, fragment_data in fragments_data.items():
        if fragment_data["absconf_index"] != 0:
            continue
        molgraph.set_nodes_visible(visible_nodes=fragment_data["visible_nodes"])
        geometry = fragment_data["final_geometry"]
        global_indices = np.array(fragment_data["visible_nodes"])
        prios = priorities[global_indices]
        core_indices = fragment_data["core_indices"]
        absconf = AbsoluteConfIndex(geometry, prios, core_indices)
        G0s = absconf.get_absconf_index()
        fragment_data["absconf_index"] = G0s
    molgraph.set_nodes_visible(include_all=True)
    return fragments_data

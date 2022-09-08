import time
import numpy as np
import pandas as pd
import logging

import trackml.dataset

from .odd_tools import load_hits_and_truth_as_trackml


#####################################################
#                   UTILD PANDAS                    #
#####################################################
def select_min(test_val, current_val):
    return min(test_val, current_val)


def select_max(test_val, current_val):
    if current_val == -1:
        return test_val
    else:
        return max(test_val, current_val)


def find_channel0_min(cells_in, nb_hits):
    cell_idx = cells_in.index.values.reshape(-1, 1)
    cells = cells_in[["hit_id", "channel0"]].values
    where_min = find_channel0_property(cells, nb_hits, select_min, 10**8)
    return where_min


def find_channel0_max(cells_in, nb_hits):
    cells = cells_in[["hit_id", "channel0"]].values
    where_max = find_channel0_property(cells, nb_hits, select_max, -(10**8))
    return where_max


def find_channel0_property(cells, nb_hits, comparator, init_val):
    nb_cells = cells.shape[0]
    cells = sort_cells_by_hit_id(cells)

    hit_property = [init_val] * nb_hits
    cell_property = [0] * nb_cells
    cell_values = cells[:, 2].tolist()
    hit_ids = cells[:, 1].tolist()

    hit_property_id = 0
    current_hit_id = hit_ids[0]
    for i, (h, v) in enumerate(zip(hit_ids, cell_values)):
        if h > current_hit_id:
            hit_property_id += 1
            current_hit_id = h
        hit_property[hit_property_id] = comparator(v, hit_property[hit_property_id])

    hit_property_id = 0
    current_hit_id = hit_ids[0]
    for i, (h, v) in enumerate(zip(hit_ids, cell_values)):
        if h > current_hit_id:
            hit_property_id += 1
            current_hit_id = h
        if v == hit_property[hit_property_id]:
            cell_property[i] = 1

    original_order = np.argsort(cells[:, 0])
    cell_property = np.array(cell_property, dtype=bool)[original_order]
    return cell_property


def sort_cells_by_hit_id(cells):
    orig_order = np.arange(len(cells)).reshape(-1, 1)
    cells = np.concatenate((orig_order, cells), 1)
    sort_idx = np.argsort(cells[:, 1])  # Sort by hit ID
    cells = cells[sort_idx]
    return cells


#################################################
#                   EXTRACT DIR                 #
#################################################


def local_angle(cell, module):
    n_u = max(cell["channel0"]) - min(cell["channel0"]) + 1
    n_v = max(cell["channel1"]) - min(cell["channel1"]) + 1
    l_u = n_u * module.pitch_u.values  # x
    l_v = n_v * module.pitch_v.values  # y
    l_w = 2 * module.module_t.values  # z
    return (l_u, l_v, l_w)


def extract_rotation_matrix(module):
    rot_matrix = np.matrix(
        [
            [module.rot_xu.values[0], module.rot_xv.values[0], module.rot_xw.values[0]],
            [module.rot_yu.values[0], module.rot_yv.values[0], module.rot_yw.values[0]],
            [module.rot_zu.values[0], module.rot_zv.values[0], module.rot_zw.values[0]],
        ]
    )
    return rot_matrix, np.linalg.inv(rot_matrix)


def cartesion_to_spherical(x, y, z):
    r3 = np.sqrt(x**2 + y**2 + z**2)
    phi = np.arctan2(y, x)
    theta = np.arccos(z / r3)
    return r3, theta, phi


def theta_to_eta(theta):
    return -np.log(np.tan(0.5 * theta))


def get_all_local_angles(hits, cells, detector):
    direction_count_u = cells.groupby(["hit_id"]).channel0.agg(["min", "max"])
    direction_count_v = cells.groupby(["hit_id"]).channel1.agg(["min", "max"])
    nb_u = direction_count_u["max"] - direction_count_u["min"] + 1
    nb_v = direction_count_v["max"] - direction_count_v["min"] + 1

    vols = hits["volume_id"].values
    layers = hits["layer_id"].values
    modules = hits["module_id"].values

    pitch = detector["pixel_size"]
    thickness = detector["thicknesses"]

    pitch_cells = pitch[vols, layers, modules]
    thickness_cells = thickness[vols, layers, modules]

    l_u = nb_u * pitch_cells[:, 0]
    l_v = nb_v * pitch_cells[:, 1]
    l_w = 2 * thickness_cells
    return l_u, l_v, l_w


def get_all_rotated(hits, detector, l_u, l_v, l_w):
    vols = hits["volume_id"].values
    layers = hits["layer_id"].values
    modules = hits["module_id"].values
    rotations = detector["rotations"]
    rotations_hits = rotations[vols, layers, modules]
    u = l_u.values.reshape(-1, 1)
    v = l_v.values.reshape(-1, 1)
    w = l_w.reshape(-1, 1)
    dirs = np.concatenate((u, v, w), axis=1)

    dirs = np.expand_dims(dirs, axis=2)
    vecRot = np.matmul(rotations_hits, dirs).squeeze(2)
    return vecRot


def extract_dir_new(hits, cells, detector):
    l_u, l_v, l_w = get_all_local_angles(hits, cells, detector)
    g_matrix_all = get_all_rotated(hits, detector, l_u, l_v, l_w)
    hit_ids, cell_counts, cell_vals = (
        hits["hit_id"].to_numpy(),
        hits["cell_count"].to_numpy(),
        hits["cell_val"].to_numpy(),
    )

    l_u, l_v = l_u.to_numpy(), l_v.to_numpy()

    _, g_theta, g_phi = np.vstack(cartesion_to_spherical(*list(g_matrix_all.T)))
    logging.info("G calc")
    _, l_theta, l_phi = cartesion_to_spherical(l_u, l_v, l_w)
    logging.info("L calc")
    l_eta = theta_to_eta(l_theta)
    g_eta = theta_to_eta(g_theta)

    angles = np.vstack(
        [hit_ids, cell_counts, cell_vals, l_eta, l_phi, l_u, l_v, l_w, g_eta, g_phi]
    ).T
    logging.info("Concated")
    df_angles = pd.DataFrame(
        angles,
        columns=[
            "hit_id",
            "cell_count",
            "cell_val",
            "leta",
            "lphi",
            "lx",
            "ly",
            "lz",
            "geta",
            "gphi",
        ],
    )
    logging.info("DF constructed")

    return df_angles


def check_diff(h1, h2, name):
    n1 = h1[name].values
    n2 = h2[name].values
    print(name, max(np.absolute(n1 - n2)))


#############################################
#           FEATURE_AUGMENTATION            #
#############################################


def augment_hit_features(hits, cells, detector_orig, detector_proc):

    cell_stats = get_cell_stats(cells)
    hits["cell_count"] = cell_stats[:, 0]
    hits["cell_val"] = cell_stats[:, 1]

    angles = extract_dir_new(hits, cells, detector_proc)

    return angles


def get_cell_stats(cells):
    hit_cells = cells.groupby(["hit_id"]).value.count().values
    hit_value = cells.groupby(["hit_id"]).value.sum().values
    cell_stats = np.hstack((hit_cells.reshape(-1, 1), hit_value.reshape(-1, 1)))
    cell_stats = cell_stats.astype(np.float32)
    return cell_stats


###########################################
#               EVENT LOADING             #
###########################################


def get_one_event(event_path, detector_orig, detector_proc):

    logging.info("Loading ODD data")
    
    hits, _ = load_hits_and_truth_as_trackml(event_path, detector_orig, mask_simhits=True)
    cells = pd.read_csv(event_path + "-cells.csv")

    logging.info("Hits and cells retrieved")

    if True:
    #try:
        angles = augment_hit_features(hits, cells, detector_orig, detector_proc)
    #except Exception as e:
    #    print(e)
    #    raise Exception("Error augmenting hits.")

    return angles

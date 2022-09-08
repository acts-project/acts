"""Utilities for processing the overall event.

The module contains useful functions for handling data at the event level. More fine-grained utilities are 
reserved for `detector_utils` and `cell_utils`.
    
Todo:
    * Pull module IDs out into a csv file for readability """

# System
import os
import argparse
import logging
import multiprocessing as mp
from functools import partial

# Externals
import yaml
import numpy as np
import pandas as pd
import trackml.dataset

import torch
from torch_geometric.data import Data

from itertools import permutations
import itertools

# Locals
from .cell_utils import get_one_event
from .odd_tools import load_hits_and_truth_as_trackml


def get_cell_information(
    data, cell_features, detector_orig, detector_proc, endcaps, noise
):

    event_file = data.event_file
    evtid = event_file[-4:]

    angles = get_one_event(event_file, detector_orig, detector_proc)
    logging.info("Angles: {}".format(angles))
    hid = pd.DataFrame(data.hid.numpy(), columns=["hit_id"])
    cell_data = torch.from_numpy(
        (hid.merge(angles, on="hit_id")[cell_features]).to_numpy()
    ).float()
    logging.info("DF merged")
    data.cell_data = cell_data

    return data


def get_layerwise_edges(hits):

    hits = hits.assign(
        R=np.sqrt(
            (hits.x - hits.vx) ** 2 + (hits.y - hits.vy) ** 2 + (hits.z - hits.vz) ** 2
        )
    )
    hits = hits.sort_values("R").reset_index(drop=True).reset_index(drop=False)
    hits.loc[hits["particle_id"] == 0, "particle_id"] = np.nan
    hit_list = (
        hits.groupby(["particle_id", "layer"], sort=False)["index"]
        .agg(lambda x: list(x))
        .groupby(level=0)
        .agg(lambda x: list(x))
    )

    true_edges = []
    for row in hit_list.values:
        for i, j in zip(row[0:-1], row[1:]):
            true_edges.extend(list(itertools.product(i, j)))
    true_edges = np.array(true_edges).T

    return true_edges, hits


def get_modulewise_edges(hits):

    signal = hits[
        ((~hits.particle_id.isna()) & (hits.particle_id != 0)) & (~hits.vx.isna())
    ]
    signal = signal.drop_duplicates(
        subset=["particle_id", "volume_id", "layer_id", "module_id"]
    )

    # Sort by increasing distance from production
    signal = signal.assign(
        R=np.sqrt(
            (signal.x - signal.vx) ** 2
            + (signal.y - signal.vy) ** 2
            + (signal.z - signal.vz) ** 2
        )
    )
    signal = signal.sort_values("R").reset_index(drop=False)

    # Handle re-indexing
    signal = signal.rename(columns={"index": "unsorted_index"}).reset_index(drop=False)
    signal.loc[signal["particle_id"] == 0, "particle_id"] = np.nan

    # Group by particle ID
    signal_list = signal.groupby(["particle_id"], sort=False)["index"].agg(
        lambda x: list(x)
    )

    true_edges = []
    for row in signal_list.values:
        for i, j in zip(row[:-1], row[1:]):
            true_edges.append([i, j])

    true_edges = np.array(true_edges).T

    true_edges = signal.unsorted_index.values[true_edges]

    return true_edges


def select_hits(hits, truth, particles, endcaps=False, noise=False):
    # TODO add mechanism to insert ignore volume-layer-ids
    if not endcaps:
        raise RuntimeError("WARNING: endcaps option not supported currently")

    # Select barrel layers and assign convenient layer number [0-9]
    vlid_groups = hits.groupby(["volume_id", "layer_id"])
    vlids = list(vlid_groups.groups.keys())
    n_det_layers = len(vlids)

    hits = pd.concat(
        [vlid_groups.get_group(vlids[i]).assign(layer=i) for i in range(n_det_layers)]
    )

    if noise:
        truth = truth.merge(
            particles[["particle_id", "vx", "vy", "vz"]], on="particle_id", how="left"
        )
    else:
        truth = truth.merge(
            particles[["particle_id", "vx", "vy", "vz"]], on="particle_id", how="inner"
        )

    truth = truth.assign(pt=np.sqrt(truth.tpx**2 + truth.tpy**2))

    # Calculate derived hits variables
    r = np.sqrt(hits.x**2 + hits.y**2)
    phi = np.arctan2(hits.y, hits.x)
    # Select the data columns we need
    hits = hits.assign(r=r, phi=phi).merge(truth, on="hit_id")

    return hits


def build_event(
    event_file,
    feature_scale,
    endcaps=False,
    modulewise=True,
    layerwise=True,
    noise=False,
    detector=None,
    cell_information=False,
    truth_hits=True,
):
    # Import basic information
    particles = pd.read_csv(event_file + "-particles.csv")
    particles["pt"] = np.sqrt(particles.px**2 + particles.py**2)
    
    truth = pd.read_csv(event_file + "-truth.csv")
    truth = truth.rename(columns={"tx": "x", "ty": "y", "tz": "z", "tt": "t"})
    truth = truth.merge(
        detector[["geometry_id", "volume_id", "layer_id", "module_id"]], on="geometry_id"
    )
    truth = truth.merge(
        particles[["particle_id", "vx", "vy", "vz", "pt"]], on="particle_id"
    )
    truth["weight"] = np.ones(len(truth.index)) / len(truth.index)
    
    # Select wether to use true hits or measurements
    if truth_hits:
        logging.info("Using truth hit information")
        hits = truth
        hits["hit_id"] = np.array(hits.index)+1
        hits = hits.drop("index", 1)
    else:     
        logging.info("Using measurement information")
        measurements = pd.read_csv(event_file + "-measurements.csv")
        simhit_map = pd.read_csv(event_file + "-measurement-simhit-map.csv")
        
        global_pos = local_to_global(measurements, detector)
        
        hits = pd.DataFrame()
        hits["geometry_id"] = measurements.geometry_id
        hits["x"] = global_pos[:,0]
        hits["y"] = global_pos[:,1]
        hits["z"] = global_pos[:,2]
        hits["hit_id"] = measurements["measurement_id"].map(dict(zip(simhit_map.measurement_id, simhit_map.hit_id)))
        hits["particle_id"] = hits["hit_id"].map(dict(zip(truth.index, truth.particle_id)))
    
    # Compute cylinder coordinates
    hits["r"] = np.sqrt(hits.x**2 + hits.y**2)
    hits["phi"] = np.arctan2(hits.y, hits.x)
    
    # Make a unique module ID and attach to hits
    if detector is not None:
        module_lookup = detector.reset_index()[["index", "volume_id", "layer_id", "module_id"]].rename(columns={"index": "module_index"})
        hits = hits.merge(module_lookup, on=["volume_id", "layer_id", "module_id"], how="left")
        module_id = hits.module_index.to_numpy()
    else:
        module_id = None

    # Seems not to be needed right now
    try:
        layer_id = hits.layer.to_numpy()
    except:
        layer_id = None
    
    # Seemst not to be needed right now / anymore
    if False:
        hits = select_hits(
            hits, truth, particles, endcaps=endcaps, noise=noise
        ).assign(evtid=int(event_file[-9:]))

    # Handle which truth graph(s) are being produced
    modulewise_true_edges, layerwise_true_edges = None, None

    if layerwise:
        layerwise_true_edges, hits = get_layerwise_edges(hits)
        logging.info(
            "Layerwise truth graph built for {} with size {}".format(
                event_file, layerwise_true_edges.shape
            )
        )
            
        pid = hits.particle_id.to_numpy()  
        assert (pid[layerwise_true_edges[0,:]] == pid[layerwise_true_edges[1,:]]).all()

    if modulewise:
        modulewise_true_edges = get_modulewise_edges(hits)
        logging.info(
            "Modulewise truth graph built for {} with size {}".format(
                event_file, modulewise_true_edges.shape
            )
        )
            
        pid = hits.particle_id.to_numpy()
        assert (pid[modulewise_true_edges[0,:]] == pid[modulewise_true_edges[1,:]]).all()

    edge_weights = (
        hits.weight.to_numpy()[modulewise_true_edges]
        if modulewise
        else hits.weight.to_numpy()[layerwise_true_edges]
    )
    edge_weight_average = (edge_weights[0] + edge_weights[1]) / 2
    edge_weight_norm = edge_weight_average / edge_weight_average.mean()

    logging.info("Weights constructed")

    return (
        hits[["r", "phi", "z"]].to_numpy() / feature_scale,
        hits.particle_id.to_numpy(),
        layer_id,
        module_id,
        modulewise_true_edges,
        layerwise_true_edges,
        hits["hit_id"].to_numpy(),
        hits.pt.to_numpy(),
        edge_weight_norm,
    )


def prepare_event(
    event_file,
    detector_orig,
    detector_proc,
    cell_features,
    progressbar=None,
    output_dir=None,
    endcaps=False,
    modulewise=True,
    layerwise=True,
    noise=False,
    cell_information=True,
    overwrite=False,
    **kwargs
):
    if True:
    #try:
        evtid = int(event_file[-9:])
        filename = os.path.join(output_dir, str(evtid))

        if not os.path.exists(filename) or overwrite:
            logging.info("Preparing event {}".format(evtid))
            feature_scale = [1000, np.pi, 1000]

            (
                X,
                pid,
                layer_id,
                module_id,
                modulewise_true_edges,
                layerwise_true_edges,
                hid,
                pt,
                weights,
            ) = build_event(
                event_file,
                feature_scale,
                endcaps=endcaps,
                modulewise=modulewise,
                layerwise=layerwise,
                noise=noise,
                detector=detector_orig,
                cell_information=cell_information
            )

            data = Data(
                x=torch.from_numpy(X).float(),
                pid=torch.from_numpy(pid),
                modules=torch.from_numpy(module_id),
                event_file=event_file,
                hid=torch.from_numpy(hid),
                pt=torch.from_numpy(pt),
                weights=torch.from_numpy(weights),
            )
            if modulewise_true_edges is not None:
                data.modulewise_true_edges = torch.from_numpy(modulewise_true_edges)
            if layerwise_true_edges is not None:
                data.layerwise_true_edges = torch.from_numpy(layerwise_true_edges)
            logging.info("Getting cell info")

            if cell_information:
                data = get_cell_information(
                    data, cell_features, detector_orig, detector_proc, endcaps, noise
                )

            with open(filename, "wb") as pickle_file:
                torch.save(data, pickle_file)

        else:
            logging.info("{} already exists".format(evtid))
    #except Exception as inst:
    #    print("File:", event_file, "had exception", inst)

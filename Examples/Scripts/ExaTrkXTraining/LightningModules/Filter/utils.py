import os
import logging
import random
from time import time as tt

import torch
import torch.nn as nn
import scipy as sp
import numpy as np
from tqdm import tqdm

from torch_geometric.data import Dataset

device = "cuda" if torch.cuda.is_available() else "cpu"

def load_dataset(input_subdir, num, **kwargs):
    if input_subdir is not None:
        all_events = os.listdir(input_subdir)
        if "sorted_events" in kwargs.keys() and kwargs["sorted_events"]:
            all_events = sorted(all_events)
        else:
            random.shuffle(all_events)
        
        all_events = [os.path.join(input_subdir, event) for event in all_events]   
        
        loaded_events = []
        for event in tqdm(all_events[:num]):
            try:
                loaded_event = torch.load(event, map_location=torch.device("cpu"))
                loaded_events.append(loaded_event)
                logging.info("Loaded event: {}".format(loaded_event.event_file))
            except:
                logging.info("Corrupted event file: {}".format(event))
        return loaded_events
    else:
        return None


def graph_intersection(pred_graph, truth_graph):
    array_size = max(pred_graph.max().item(), truth_graph.max().item()) + 1

    l1 = pred_graph.cpu().numpy()
    l2 = truth_graph.cpu().numpy()
    e_1 = sp.sparse.coo_matrix(
        (np.ones(l1.shape[1]), l1), shape=(array_size, array_size)
    ).tocsr()
    e_2 = sp.sparse.coo_matrix(
        (np.ones(l2.shape[1]), l2), shape=(array_size, array_size)
    ).tocsr()
    e_intersection = (e_1.multiply(e_2) - ((e_1 - e_2) > 0)).tocoo()

    new_pred_graph = (
        torch.from_numpy(np.vstack([e_intersection.row, e_intersection.col]))
        .long()
        .to(device)
    )
    y = e_intersection.data > 0

    return new_pred_graph, y

def sample_fake_ratio(data, fake_ratio):
    # Sample random subset of fake edges in the graph
    num_fake = (~data.y.bool()).sum()
    true_indices, fake_indices = torch.where(data.y.bool())[0], torch.where(~data.y.bool())[0]
    fake_indices = fake_indices[torch.randperm(len(fake_indices))[:int(num_fake * fake_ratio)]]
    all_indices = torch.cat([true_indices, fake_indices])
    
    return all_indices

class LargeDataset(Dataset):
    def __init__(self, root, num_events, data_name, hparams, transform=None, pre_transform=None, pre_filter=None):
        super().__init__(root, transform, pre_transform, pre_filter)
        self.hparams = hparams
        self.data_name = data_name
        
        self.input_paths = os.listdir(root)
        if "sorted_events" in hparams.keys() and hparams["sorted_events"]:
            self.input_paths = sorted(self.input_paths)
        else:
            random.shuffle(self.input_paths)
        
        self.input_paths = [os.path.join(root, event) for event in self.input_paths][:num_events]
        
    def len(self):
        return len(self.input_paths)

    def get(self, idx):
        data = torch.load(self.input_paths[idx], map_location=torch.device("cpu"))
        
        # Order edges by increasing module ID
        if "volume_id" in data.keys:
            edges_to_be_flipped = data.volume_id[data.edge_index[0]] > data.volume_id[data.edge_index[1]]
            data.edge_index[:, edges_to_be_flipped] =  data.edge_index[:, edges_to_be_flipped].flip(0)
            assert (data.volume_id[data.edge_index[0]] <= data.volume_id[data.edge_index[1]]).all(), "Flip didn't work!"
        

        if "train_fake_sample" in self.hparams.keys() and self.hparams["train_fake_sample"] and (self.data_name=="train"):
            sampled_indices = sample_fake_ratio(data, self.hparams["train_fake_sample"])
            data.edge_index = data.edge_index[:, sampled_indices]
            data.y = data.y[sampled_indices]

        data.y_pid = (data.pid[data.edge_index[0]] == data.pid[data.edge_index[1]]) & data.pid[data.edge_index[0]].bool()
        return data
    

class filter_dataset(Dataset):
    def __init__(self, dataset, hparams):

        # Setup here
        self.dataset = dataset
        self.hparams = hparams

    def __len__(self):

        return len(self.dataset)

    def __getitem__(self, idx):

        batch = self.dataset[idx]

        if "subset" in self.hparams["regime"]:
            subset_mask = np.isin(
                batch.edge_index, batch.signal_true_edges.unique()
            ).any(0)
            batch.edge_index = batch.edge_index[:, subset_mask]
            batch.y = batch.y[subset_mask]

        if self.hparams["ratio"] != 0:
            num_true, num_false = (
                batch.signal_true_edges.shape[1],
                (~batch.y.bool()).sum(),
            )

            # Select a subset of fake edges randomly
            start_index = torch.randint(
                len(batch.y) - 2 * self.hparams["ratio"] * num_true, (1,)
            )
            end_index = start_index + 2 * self.hparams["ratio"] * num_true
            random_edges = batch.edge_index[:, start_index:end_index]
            combined_edges = torch.cat(
                [
                    batch.signal_true_edges,
                    batch.signal_true_edges.flip(0),
                    random_edges,
                ],
                dim=1,
            )
            combined_y = torch.cat(
                [
                    torch.ones(2 * batch.signal_true_edges.shape[1]),
                    batch.y[start_index:end_index],
                ]
            )

            # Shuffle in true edges
            shuffle_indices = torch.randperm(combined_edges.shape[1])
            combined_edges = combined_edges[:, shuffle_indices]
            combined_y = combined_y[shuffle_indices]

            # Select a further subset in order to handle memory issues
            start_index = torch.randint(
                len(combined_y) - self.hparams["edges_per_batch"], (1,)
            )
            end_index = start_index + self.hparams["edges_per_batch"]

            combined_edges = combined_edges[:, start_index:end_index]
            combined_y = combined_y[start_index:end_index]

        subbatch = {
            "x": batch.x,
            "cell_data": batch.cell_data,
            "edge_index": combined_edges,
            "y": combined_y,
        }

        return subbatch

def make_mlp(
    input_size,
    sizes,
    hidden_activation="ReLU",
    output_activation="ReLU",
    layer_norm=False,
    batch_norm=False,
):
    """Construct an MLP with specified fully-connected layers."""
    hidden_activation = getattr(nn, hidden_activation)
    if output_activation is not None:
        output_activation = getattr(nn, output_activation)
    layers = []
    n_layers = len(sizes)
    sizes = [input_size] + sizes
    # Hidden layers
    for i in range(n_layers - 1):
        layers.append(nn.Linear(sizes[i], sizes[i + 1]))
        if layer_norm:
            layers.append(nn.LayerNorm(sizes[i + 1], elementwise_affine=False))
        if batch_norm:
            layers.append(nn.BatchNorm1d(sizes[i + 1], track_running_stats=False, affine=False))
        layers.append(hidden_activation())
    # Final layer
    layers.append(nn.Linear(sizes[-2], sizes[-1]))
    if output_activation is not None:
        if layer_norm:
            layers.append(nn.LayerNorm(sizes[-1], elementwise_affine=False))
        if batch_norm:
            layers.append(nn.BatchNorm1d(sizes[-1], track_running_stats=False, affine=False))
        layers.append(output_activation())
    return nn.Sequential(*layers)

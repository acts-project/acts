# System imports
import sys
import os
import copy

# 3rd party imports
import pytorch_lightning as pl
from pytorch_lightning.callbacks import Callback
from torch.nn import Linear
import torch.nn as nn
import torch.nn.functional as F
from torch_cluster import radius_graph
import torch
from torch_geometric.data import DataLoader

# Local imports
from ..utils import graph_intersection
from ..filter_base import FilterBase, FilterBaseBalanced


class VanillaFilter(FilterBaseBalanced):
    def __init__(self, hparams):
        super().__init__(hparams)
        """
        Initialise the Lightning Module that can scan over different filter training regimes
        """

        # Construct the MLP architecture
        self.input_layer = Linear(
            hparams["in_channels"] * 2 + hparams["emb_channels"] * 2, hparams["hidden"]
        )
        layers = [
            Linear(hparams["hidden"], hparams["hidden"])
            for _ in range(hparams["nb_layer"] - 1)
        ]
        self.layers = nn.ModuleList(layers)
        self.output_layer = nn.Linear(hparams["hidden"], 1)
        self.layernorm = nn.LayerNorm(hparams["hidden"])
        self.batchnorm = nn.BatchNorm1d(
            num_features=hparams["hidden"], track_running_stats=False
        )
        self.act = nn.Tanh()

    def forward(self, x, e, emb=None):
        if emb is not None:
            x = self.input_layer(
                torch.cat([x[e[0]], emb[e[0]], x[e[1]], emb[e[1]]], dim=-1)
            )
        else:
            #x = self.input_layer(torch.cat([x[e[0]], x[e[1]]], dim=-1))
            x = self.input_layer(torch.cat(
                [
                    torch.index_select(x, 0, torch.select(e, 0, 0)),
                    torch.index_select(x, 0, torch.select(e, 0, 1)),
                ], dim=-1))

        for l in self.layers:
            x = l(x)
            x = self.act(x)
            if self.hparams["layernorm"]:
                x = self.layernorm(x)  # Option of LayerNorm
            if self.hparams["batchnorm"]:
                x = self.batchnorm(x)  # Option of Batch
        x = self.output_layer(x)
        return x


class FilterInferenceCallback(Callback):
    def __init__(self):
        self.output_dir = None
        self.overwrite = False

    def on_train_start(self, trainer, pl_module):
        # Prep the directory to produce inference data to
        self.output_dir = pl_module.hparams.output_dir
        self.datatypes = ["train", "val", "test"]
        os.makedirs(self.output_dir, exist_ok=True)
        [
            os.makedirs(os.path.join(self.output_dir, datatype), exist_ok=True)
            for datatype in self.datatypes
        ]

    def on_train_end(self, trainer, pl_module):
        """
        This method shouldn't need to change between stages
        """

        print("Training finished, running inference to filter graphs...")

        # By default, the set of examples propagated through the pipeline will be train+val+test set
        datasets = {
            "train": pl_module.trainset,
            "val": pl_module.valset,
            "test": pl_module.testset,
        }
        total_length = sum([len(dataset) for dataset in datasets.values()])
        batch_incr = 0

        pl_module.eval()
        with torch.no_grad():
            for set_idx, (datatype, dataset) in enumerate(datasets.items()):
                for batch_idx, batch in enumerate(dataset):
                    percent = (batch_incr / total_length) * 100
                    sys.stdout.flush()
                    sys.stdout.write(f"{percent:.01f}% inference complete \r")
                    if (
                        not os.path.exists(
                            os.path.join(
                                self.output_dir, datatype, batch.event_file[-4:]
                            )
                        )
                    ) or self.overwrite:
                        batch_to_save = copy.deepcopy(batch)
                        batch_to_save = batch_to_save.to(
                            pl_module.device
                        )  # Is this step necessary??
                        batch_to_save = self.construct_downstream(
                            batch_to_save, pl_module
                        ).to("cpu")
                        self.save_downstream(batch_to_save, pl_module, datatype)

                    batch_incr += 1

    def construct_downstream(self, batch, pl_module):

        """
        This contains the bulk of pipeline logic for this stage
        """
        emb = (
            None if (pl_module.hparams["emb_channels"] == 0) else batch.embedding
        )  # Does this work??

        sections = 8
        cut_list = []
        for j in range(sections):
            subset_ind = torch.chunk(torch.arange(batch.edge_index.shape[1]), sections)[
                j
            ]
            output = (
                pl_module(
                    torch.cat([batch.cell_data, batch.x], axis=-1),
                    batch.edge_index[:, subset_ind],
                    emb,
                ).squeeze()
                if ("ci" in pl_module.hparams["regime"])
                else pl_module(batch.x, batch.edge_index[:, subset_ind], emb).squeeze()
            )
            cut = F.sigmoid(output) > pl_module.hparams["filter_cut"]
            cut_list.append(cut)

        cut_list = torch.cat(cut_list)

        if "pid" not in pl_module.hparams["regime"]:
            batch.y = batch.y[cut_list]

        y_pid = batch.pid[batch.edge_index[0]] == batch.pid[batch.edge_index[1]]
        batch.y_pid = y_pid[cut_list]
        batch.edge_index = batch.edge_index[:, cut_list]
        if "weighting" in pl_module.hparams["regime"]:
            batch.weights = batch.weights[cut_list]

        return batch

    def save_downstream(self, batch, pl_module, datatype):

        with open(
            os.path.join(self.output_dir, datatype, batch.event_file[-4:]), "wb"
        ) as pickle_file:
            torch.save(batch, pickle_file)

import sys

import pytorch_lightning as pl
from pytorch_lightning import LightningModule
from pytorch_lightning.callbacks import Callback
import torch.nn as nn
from torch.nn import Linear
import torch.nn.functional as F
import torch
from torch_scatter import scatter_add, scatter_mean, scatter_max
from torch.utils.checkpoint import checkpoint

from ..gnn_base import GNNBase
from ..utils import make_mlp


class VanillaGCN(GNNBase):
    def __init__(self, hparams):
        super().__init__(hparams)
        """
        Initialise the Lightning Module that can scan over different GNN training regimes
        """

        # Setup input network
        self.node_encoder = make_mlp(
            hparams["in_channels"],
            [hparams["hidden"]] * hparams["nb_node_layer"],
            output_activation=hparams["hidden_activation"],
            layer_norm=hparams["layernorm"],
        )

        # The edge network computes new edge features from connected nodes
        self.edge_network = make_mlp(
            2 * (hparams["hidden"]),
            [hparams["hidden"]] * hparams["nb_edge_layer"] + [1],
            layer_norm=hparams["layernorm"],
            output_activation=None,
            hidden_activation=hparams["hidden_activation"],
        )

        # The node network computes new node features
        self.node_network = make_mlp(
            (hparams["hidden"]) * 2,
            [hparams["hidden"]] * hparams["nb_node_layer"],
            layer_norm=hparams["layernorm"],
            hidden_activation=hparams["hidden_activation"],
        )

    def forward(self, x, edge_index):

        input_x = x

        x = self.node_encoder(x)
        #         x = F.softmax(x, dim=-1)

        start, end = edge_index

        for i in range(self.hparams["n_graph_iters"]):

            #             x_initial = x

            messages = scatter_add(
                x[start], end, dim=0, dim_size=x.shape[0]
            ) + scatter_add(x[end], start, dim=0, dim_size=x.shape[0])

            node_inputs = torch.cat([x, messages], dim=-1)
            #             node_inputs = F.softmax(node_inputs, dim=-1)

            x = self.node_network(node_inputs)

        #             x = x + x_initial

        edge_inputs = torch.cat([x[start], x[end]], dim=1)
        return self.edge_network(edge_inputs)

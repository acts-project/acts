import sys

import pytorch_lightning as pl
from pytorch_lightning import LightningModule
from pytorch_lightning.callbacks import Callback
import torch.nn as nn
from torch.nn import Linear
import torch.nn.functional as F
import torch
from torch_scatter import scatter_add
from torch.utils.checkpoint import checkpoint

from ..gnn_base import GNNBase
from ..utils import make_mlp


class VanillaCheckResAGNN(GNNBase):
    def __init__(self, hparams):
        super().__init__(hparams)
        """
        Initialise the Lightning Module that can scan over different GNN training regimes
        """
        self.edge_network = make_mlp(
            (hparams["in_channels"] + hparams["hidden"]) * 2,
            [hparams["in_channels"] + hparams["hidden"]] * hparams["nb_edge_layer"]
            + [1],
            hidden_activation=hparams["hidden_activation"],
            output_activation=None,
            layer_norm=hparams["layernorm"],
        )

        self.node_network = make_mlp(
            (hparams["in_channels"] + hparams["hidden"]) * 2,
            [hparams["hidden"]] * hparams["nb_node_layer"],
            hidden_activation=hparams["hidden_activation"],
            output_activation=None,
            layer_norm=hparams["layernorm"],
        )

        self.input_network = make_mlp(
            hparams["in_channels"],
            [hparams["hidden"]] * hparams["nb_node_layer"],
            output_activation=hparams["hidden_activation"],
            layer_norm=hparams["layernorm"],
        )

    def forward(self, x, edge_index):
        start, end = edge_index
        input_x = x

        x = self.input_network(x)

        # Shortcut connect the inputs onto the hidden representation
        x = torch.cat([x, input_x], dim=-1)

        # Loop over iterations of edge and node networks
        for i in range(self.hparams["n_graph_iters"]):
            x_inital = x

            # Apply edge network
            edge_inputs = torch.cat([x[start], x[end]], dim=1)
            e = checkpoint(self.edge_network, edge_inputs)
            e = torch.sigmoid(e)

            # Apply node network
            messages = scatter_add(
                e * x[start], end, dim=0, dim_size=x.shape[0]
            ) + scatter_add(e * x[end], start, dim=0, dim_size=x.shape[0])
            node_inputs = torch.cat([messages, x], dim=1)
            x = checkpoint(self.node_network, node_inputs)

            # Shortcut connect the inputs onto the hidden representation
            x = torch.cat([x, input_x], dim=-1)

            # Residual connection
            x = x_inital + x

        edge_inputs = torch.cat([x[start], x[end]], dim=1)
        return checkpoint(self.edge_network, edge_inputs)

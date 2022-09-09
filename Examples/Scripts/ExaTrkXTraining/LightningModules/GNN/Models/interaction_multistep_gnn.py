import sys

import pytorch_lightning as pl
from pytorch_lightning import LightningModule
from pytorch_lightning.callbacks import Callback
import torch.nn as nn
from torch.nn import Linear
import torch.nn.functional as F
import torch
from torch_scatter import scatter_add
from torch_geometric.nn.conv import MessagePassing
from torch.utils.checkpoint import checkpoint

from ..gnn_base import GNNBase
from ..utils import make_mlp


class InteractionMultistepGNN(GNNBase):
    def __init__(self, hparams):
        super().__init__(hparams)
        """
        Initialise the Lightning Module that can scan over different GNN training regimes
        """

        # Setup input network
        self.node_encoder = make_mlp(
            hparams["in_channels"],
            [hparams["hidden"]],
            output_activation=hparams["hidden_activation"],
            layer_norm=hparams["layernorm"],
        )

        # The edge network computes new edge features from connected nodes
        self.edge_encoder = make_mlp(
            2 * (hparams["hidden"]),
            [hparams["hidden"]] * hparams["nb_edge_layer"],
            layer_norm=hparams["layernorm"],
            output_activation=None,
            hidden_activation=hparams["hidden_activation"],
        )

        # The edge network computes new edge features from connected nodes
        self.edge_network = make_mlp(
            4 * hparams["hidden"],
            [hparams["hidden"]] * hparams["nb_edge_layer"],
            layer_norm=hparams["layernorm"],
            output_activation=None,
            hidden_activation=hparams["hidden_activation"],
        )

        # The node network computes new node features
        self.node_network = make_mlp(
            4 * hparams["hidden"],
            [hparams["hidden"]] * hparams["nb_node_layer"],
            layer_norm=hparams["layernorm"],
            output_activation=None,
            hidden_activation=hparams["hidden_activation"],
        )

        # Final edge output classification network
        self.output_edge_classifier = make_mlp(
            3 * hparams["hidden"],
            [hparams["hidden"], 1],
            layer_norm=hparams["layernorm"],
            output_activation=None,
            hidden_activation=hparams["hidden_activation"],
        )

    def forward(self, x, edge_index):

        start, end = edge_index

        # Encode the graph features into the hidden space
        x = self.node_encoder(x)
        e = self.edge_encoder(torch.cat([x[start], x[end]], dim=1))
        input_x = x
        input_e = e

        edge_outputs = []
        # Loop over iterations of edge and node networks
        for i in range(self.hparams["n_graph_iters"]):

            # Cocnatenate with initial latent space
            x = torch.cat([x, input_x], dim=-1)
            e = torch.cat([e, input_e], dim=-1)

            # Compute new node features
            edge_messages = scatter_add(
                e, end, dim=0, dim_size=x.shape[0]
            ) + scatter_add(e, start, dim=0, dim_size=x.shape[0])
            node_inputs = torch.cat([x, edge_messages], dim=-1)
            x = self.node_network(node_inputs)

            # Compute new edge features
            edge_inputs = torch.cat([x[start], x[end], e], dim=-1)
            e = self.edge_network(edge_inputs)
            e = torch.sigmoid(e)

            classifier_inputs = torch.cat([x[start], x[end], e], dim=1)
            edge_outputs.append(
                self.output_edge_classifier(classifier_inputs).squeeze(-1)
            )

        # Compute final edge scores; use original edge directions only
        return edge_outputs

    def training_step(self, batch, batch_idx):

        weight = (
            torch.tensor(self.hparams["weight"])
            if ("weight" in self.hparams)
            else torch.tensor((~batch.y_pid.bool()).sum() / batch.y_pid.sum())
        )

        output = (
            self(torch.cat([batch.cell_data, batch.x], axis=-1), batch.edge_index)
            if ("ci" in self.hparams["regime"])
            else self(batch.x, batch.edge_index)
        )

        if "pid" in self.hparams["regime"]:
            y_pid = (
                batch.pid[batch.edge_index[0]] == batch.pid[batch.edge_index[1]]
            ).float()
            y_pid = y_pid.repeat((self.hparams["n_graph_iters"]))
            loss = F.binary_cross_entropy_with_logits(
                torch.cat(output), y_pid.float(), pos_weight=weight
            )
        else:
            y = batch.y.repeat((self.hparams["n_graph_iters"]))
            loss = F.binary_cross_entropy_with_logits(
                torch.cat(output), y, pos_weight=weight
            )

        result = pl.TrainResult(minimize=loss)
        result.log("train_loss", loss, prog_bar=True)

        return result

    def validation_step(self, batch, batch_idx):

        weight = (
            torch.tensor(self.hparams["weight"])
            if ("weight" in self.hparams)
            else torch.tensor((~batch.y_pid.bool()).sum() / batch.y_pid.sum())
        )

        output = (
            self(torch.cat([batch.cell_data, batch.x], axis=-1), batch.edge_index)
            if ("ci" in self.hparams["regime"])
            else self(batch.x, batch.edge_index)
        )

        if "pid" in self.hparams["regime"]:
            y_pid = (
                batch.pid[batch.edge_index[0]] == batch.pid[batch.edge_index[1]]
            ).float()
            val_loss = F.binary_cross_entropy_with_logits(
                torch.cat(output),
                y_pid.float().repeat((self.hparams["n_graph_iters"])),
                pos_weight=weight,
            )
        else:
            y = batch.y
            val_loss = F.binary_cross_entropy_with_logits(
                torch.cat(output),
                y.repeat((self.hparams["n_graph_iters"])),
                pos_weight=weight,
            )

        result = pl.EvalResult(checkpoint_on=val_loss)
        result.log("val_loss", val_loss)

        # Edge filter performance
        preds = F.sigmoid(output[-1]) > self.hparams["edge_cut"]  # Maybe send to CPU??
        edge_positive = preds.sum().float()

        if "pid" in self.hparams["regime"]:
            edge_true = y_pid.sum().float()
            edge_true_positive = (y_pid.bool() & preds).sum().float()
        else:
            edge_true = y.sum()
            edge_true_positive = (y.bool() & preds).sum().float()

        result.log_dict(
            {
                "eff": torch.tensor(edge_true_positive / edge_true),
                "pur": torch.tensor(edge_true_positive / edge_positive),
            }
        )

        return result


class CheckpointedInteractionMultistepGNN(InteractionMultistepGNN):
    def __init__(self, hparams):
        super().__init__(hparams)
        """
        Initialise the Lightning Module that can scan over different GNN training regimes
        """

    def forward(self, x, edge_index):

        start, end = edge_index

        # Encode the graph features into the hidden space
        x = self.node_encoder(x)
        e = self.edge_encoder(torch.cat([x[start], x[end]], dim=1))
        input_x = x
        input_e = e

        edge_outputs = []
        # Loop over iterations of edge and node networks
        for i in range(self.hparams["n_graph_iters"]):

            # Cocnatenate with initial latent space
            x = torch.cat([x, input_x], dim=-1)
            e = torch.cat([e, input_e], dim=-1)

            # Compute new node features
            edge_messages = scatter_add(
                e, end, dim=0, dim_size=x.shape[0]
            ) + scatter_add(e, start, dim=0, dim_size=x.shape[0])
            node_inputs = torch.cat([x, edge_messages], dim=-1)
            x = checkpoint(self.node_network, node_inputs)

            # Compute new edge features
            edge_inputs = torch.cat([x[start], x[end], e], dim=-1)
            e = checkpoint(self.edge_network, edge_inputs)
            e = torch.sigmoid(e)

            classifier_inputs = torch.cat([x[start], x[end], e], dim=1)
            edge_outputs.append(
                checkpoint(self.output_edge_classifier, classifier_inputs).squeeze(-1)
            )

        # Compute final edge scores; use original edge directions only
        return edge_outputs

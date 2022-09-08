import sys, os
import logging

import pytorch_lightning as pl
from pytorch_lightning import LightningModule
from datetime import timedelta
import torch.nn.functional as F
from torch_geometric.data import DataLoader
from torch.nn import Linear
import torch

from .utils import load_dataset, random_edge_slice_v2


class RegressionBase(LightningModule):
    def __init__(self, hparams):
        super().__init__()
        """
        Initialise the Lightning Module that can scan over different GNN training regimes
        """
        # Assign hyperparameters
        self.hparams = hparams
        self.hparams["posted_alert"] = False

    def setup(self, stage):
        # Handle any subset of [train, val, test] data split, assuming that ordering
        input_dirs = [None, None, None]
        input_dirs[: len(self.hparams["datatype_names"])] = [
            os.path.join(self.hparams["input_dir"], datatype)
            for datatype in self.hparams["datatype_names"]
        ]
        self.trainset, self.valset, self.testset = [
            load_dataset(
                input_dir, self.hparams["datatype_split"][i], self.hparams["pt_min"]
            )
            for i, input_dir in enumerate(input_dirs)
        ]

    def train_dataloader(self):
        if self.trainset is not None:
            return DataLoader(self.trainset, batch_size=1, num_workers=1)
        else:
            return None

    def val_dataloader(self):
        if self.valset is not None:
            return DataLoader(self.valset, batch_size=1, num_workers=1)
        else:
            return None

    def test_dataloader(self):
        if self.testset is not None:
            return DataLoader(self.testset, batch_size=1, num_workers=1)
        else:
            return None

    def configure_optimizers(self):
        optimizer = [
            torch.optim.AdamW(
                self.parameters(),
                lr=(self.hparams["lr"]),
                betas=(0.9, 0.999),
                eps=1e-08,
                amsgrad=True,
            )
        ]
        scheduler = [
            {
                "scheduler": torch.optim.lr_scheduler.StepLR(
                    optimizer[0],
                    step_size=self.hparams["patience"],
                    gamma=self.hparams["factor"],
                ),
                "interval": "epoch",
                "frequency": 1,
            }
        ]
        return optimizer, scheduler

    def training_step(self, batch, batch_idx):

        weight = (
            torch.tensor(self.hparams["weight"])
            if ("weight" in self.hparams)
            else torch.tensor((~batch.y_pid.bool()).sum() / batch.y_pid.sum())
        )

        node_output, _ = (
            self(torch.cat([batch.cell_data, batch.x], axis=-1), batch.edge_index)
            if ("ci" in self.hparams["regime"])
            else self(batch.x, batch.edge_index)
        )
        edge_truth = (
            (batch.pid[batch.edge_index[0]] == batch.pid[batch.edge_index[1]]).float()
            if "pid" in self.hparams["regime"]
            else batch.y
        )
        node_truth = batch.pt

        if "weighting" in self.hparams["regime"]:
            manual_weights = batch.weights
        else:
            manual_weights = None

        loss = F.mse_loss(node_output.squeeze(), node_truth.float())

        if "hybrid" in self.hparams["regime"]:
            loss += F.binary_cross_entropy_with_logits(
                output, truth.float(), weight=manual_weights, pos_weight=weight
            )

        self.log("train_loss", loss)

        return loss

    def shared_evaluation(self, batch, batch_idx):

        weight = (
            torch.tensor(self.hparams["weight"])
            if ("weight" in self.hparams)
            else torch.tensor((~batch.y_pid.bool()).sum() / batch.y_pid.sum())
        )

        node_output, edge_output = (
            self(torch.cat([batch.cell_data, batch.x], axis=-1), batch.edge_index)
            if ("ci" in self.hparams["regime"])
            else self(batch.x, batch.edge_index)
        )

        edge_truth = (
            (batch.pid[batch.edge_index[0]] == batch.pid[batch.edge_index[1]]).float()
            if "pid" in self.hparams["regime"]
            else batch.y
        )
        node_truth = batch.pt

        if "weighting" in self.hparams["regime"]:
            manual_weights = batch.weights
        else:
            manual_weights = None

        loss = F.mse_loss(node_output.squeeze(), node_truth.float())

        # Node regression performance
        # True positive defined as prediction being within 10% of true pT
        node_error = torch.abs(node_output.squeeze() - node_truth) / node_truth
        node_true_positive = (node_error < 0.1).sum()

        node_accuracy = node_true_positive / node_truth.shape[0]

        # Edge classification performance
        edge_preds = F.sigmoid(edge_output) > self.hparams["edge_cut"]

        edge_positive = edge_preds.sum().float()
        edge_true = edge_truth.sum().float()
        edge_true_positive = (edge_truth.bool() & edge_preds).sum().float()

        edge_eff = torch.tensor(edge_true_positive / edge_true)
        edge_pur = torch.tensor(edge_true_positive / edge_positive)

        current_lr = self.optimizers().param_groups[0]["lr"]

        self.log_dict(
            {
                "val_loss": loss,
                "edge_eff": edge_eff,
                "edge_pur": edge_pur,
                "node_accuracy": node_accuracy,
                "current_lr": current_lr,
            }
        )

        return {
            "loss": loss,
            "edge_preds": edge_preds.cpu().numpy(),
            "edge_truth": edge_truth.cpu().numpy(),
            "node_accuracy": node_accuracy,
        }

    def validation_step(self, batch, batch_idx):

        outputs = self.shared_evaluation(batch, batch_idx)

        return outputs["loss"]

    def test_step(self, batch, batch_idx):

        outputs = self.shared_evaluation(batch, batch_idx)

        return outputs

    def test_step_end(self, output_results):

        print("Step:", output_results)

    def test_epoch_end(self, outputs):

        print("Epoch:", outputs)

    def optimizer_step(
        self,
        epoch,
        batch_idx,
        optimizer,
        optimizer_idx,
        optimizer_closure=None,
        on_tpu=False,
        using_native_amp=False,
        using_lbfgs=False,
    ):
        # warm up lr
        if (self.hparams["warmup"] is not None) and (
            self.trainer.global_step < self.hparams["warmup"]
        ):
            lr_scale = min(
                1.0, float(self.trainer.global_step + 1) / self.hparams["warmup"]
            )
            for pg in optimizer.param_groups:
                pg["lr"] = lr_scale * self.hparams["lr"]

        # update params
        optimizer.step(closure=optimizer_closure)
        optimizer.zero_grad()

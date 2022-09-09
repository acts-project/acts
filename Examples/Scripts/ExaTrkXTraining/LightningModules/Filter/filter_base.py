# System imports
import sys, os

# 3rd party imports
import pytorch_lightning as pl
from pytorch_lightning import LightningModule
import torch
from torch.nn import Linear
import torch.nn.functional as F
from torch.utils.data import random_split
from torch.utils.data import Dataset, DataLoader
from torch_geometric.data import DataLoader as GeoLoader
import numpy as np

from sklearn.metrics import roc_auc_score

device = "cuda" if torch.cuda.is_available() else "cpu"

# Local imports
from .utils import load_dataset, filter_dataset, LargeDataset


class FilterBase(LightningModule):
    def __init__(self, hparams):
        super().__init__()
        """
        Initialise the Lightning Module that can scan over different filter training regimes
        """
        self.save_hyperparameters(hparams)

    def setup(self, stage):
        # Handle any subset of [train, val, test] data split, assuming that ordering

        input_dirs = [None, None, None]
        input_dirs[: len(self.hparams["datatype_names"])] = [
            os.path.join(self.hparams["input_dir"], datatype)
            for datatype in self.hparams["datatype_names"]
        ]
        
        if "trainset" not in self.__dict__.keys():
            self.trainset, self.valset, self.testset = [
                load_dataset(input_dir, self.hparams["datatype_split"][i], **self.hparams)
                for i, input_dir in enumerate(input_dirs)
            ]
            
        if (
            "logger" in self.__dict__.keys()
            and "_experiment" in self.logger.__dict__.keys()
        ):
            self.logger.experiment.define_metric("val_loss", summary="min")
            self.logger.experiment.define_metric("auc", summary="max")

    def train_dataloader(self):
        self.trainset = filter_dataset(self.trainset, self.hparams)
        if self.trainset is not None:
            return DataLoader(self.trainset, batch_size=1, num_workers=4)
        else:
            return None

    def val_dataloader(self):
        if self.valset is not None:
            return GeoLoader(self.valset, batch_size=1, num_workers=1)
        else:
            return None

    def test_dataloader(self):
        if self.testset is not None:
            return GeoLoader(self.testset, batch_size=1, num_workers=1)
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

    def get_input_data(self, batch):

        if "ci" in self.hparams["regime"]:
            input_data = torch.cat(
                [batch.cell_data[:, : self.hparams["cell_channels"]], batch.x], axis=-1
            )
            input_data[input_data != input_data] = 0
        else:
            input_data = batch.x
            input_data[input_data != input_data] = 0

        return input_data

    def training_step(self, batch, batch_idx):

        positive_weight = (
            torch.tensor(self.hparams["weight"])
            if ("weight" in self.hparams)
            else torch.tensor(self.hparams["ratio"])
        )

        batch["cell_data"], batch["x"], batch["edge_index"], batch["y"] = (
            batch["cell_data"].squeeze(),
            batch["x"].squeeze(),
            batch["edge_index"].squeeze(),
            batch["y"].squeeze(),
        )
        input_data = self.get_input_data(batch)

        output = self(input_data, batch["edge_index"]).squeeze()

        if "weighting" in self.hparams["regime"]:
            manual_weights = batch.weights[combined_indices]
            manual_weights[batch.y[combined_indices] == 0] = 1
        else:
            manual_weights = None

        if "pid" in self.hparams["regime"]:
            y_pid = (
                batch.pid[batch.edge_index[0, combined_indices]]
                == batch.pid[batch.edge_index[1, combined_indices]]
            )
            loss = F.binary_cross_entropy_with_logits(
                output, y_pid.float(), weight=manual_weights, pos_weight=positive_weight
            )
        else:
            loss = F.binary_cross_entropy_with_logits(
                output,
                batch["y"].float(),
                weight=manual_weights,
                pos_weight=positive_weight,
            )
        self.log("train_loss", loss)

        return loss

    def get_metrics(self, truth, scores):

        predictions = scores > self.hparams["filter_cut"]

        edge_positive = predictions.sum().float()
        edge_true = truth.sum().float()
        edge_true_positive = (truth.bool() & predictions).sum().float()

        eff = edge_true_positive / edge_true
        pur = edge_true_positive / edge_positive

        auc = roc_auc_score(truth.bool().cpu().detach(), scores.cpu().detach())

        return eff, pur, auc

    def shared_evaluation(self, batch, batch_idx, log=False):

        """
        This method is shared between validation steps and test steps
        """
        
        emb = (
            None if (self.hparams["emb_channels"] == 0) else batch.embedding
        )  # Does this work??

        score_list = []
        val_loss = torch.tensor(0).to(self.device)
        for j in range(self.hparams["n_chunks"]):
            
            subset_ind = torch.chunk(
                torch.arange(batch.edge_index.shape[1]), self.hparams["n_chunks"]
            )[j]

            if "ci" in self.hparams["regime"]:
                output = self(
                    torch.cat(
                        [batch.cell_data[:, : self.hparams["cell_channels"]], batch.x],
                        axis=-1,
                    ),
                    batch.edge_index[:, subset_ind],
                    emb,
                ).squeeze()
            else:
                output = self(batch.x, batch.edge_index[:, subset_ind], emb).squeeze()

            scores = torch.sigmoid(output)
            score_list.append(scores)

            if "weighting" in self.hparams["regime"]:
                manual_weights = batch.weights[subset_ind]
                manual_weights[batch.y[subset_ind] == 0] = 1
            else:
                manual_weights = None

            if "pid" not in self.hparams["regime"]:
                val_loss = val_loss + F.binary_cross_entropy_with_logits(
                    output, batch.y[subset_ind].float(), weight=manual_weights
                )
            else:
                y_pid = (
                    batch.pid[batch.edge_index[0, subset_ind]]
                    == batch.pid[batch.edge_index[1, subset_ind]]
                )
                val_loss = +F.binary_cross_entropy_with_logits(
                    output, y_pid.float(), weight=manual_weights
                )

        score_list = torch.cat(score_list)
        cut_list = score_list > self.hparams["edge_cut"]

        # Edge filter performance
        edge_positive = cut_list.sum().float()
        if "pid" in self.hparams["regime"]:
            true_y = batch.pid[batch.edge_index[0]] == batch.pid[batch.edge_index[1]]
        else:
            true_y = batch.y

        edge_true = true_y.sum()
        edge_true_positive = (true_y.bool() & cut_list).sum().float()

        if log:

            current_lr = self.optimizers().param_groups[0]["lr"]

            self.log_dict(
                {
                    "eff": torch.tensor(edge_true_positive / edge_true),
                    "pur": torch.tensor(edge_true_positive / edge_positive),
                    "val_loss": val_loss,
                    "current_lr": current_lr,
                }
            )
        
        return {"loss": val_loss, "preds": score_list, "truth": true_y}

    def validation_step(self, batch, batch_idx):

        outputs = self.shared_evaluation(batch, batch_idx, log=True)

        return outputs["loss"]

    def test_step(self, batch, batch_idx):
        """
        Step to evaluate the model's performance
        """
        outputs = self.shared_evaluation(batch, batch_idx, log=True)

        return outputs

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
        """
        Use this to manually enforce warm-up. In the future, this may become built-into PyLightning
        """
        # warm up lr
        if (self.hparams["warmup"] is not None) and (
            self.trainer.current_epoch < self.hparams["warmup"]
        ):
            lr_scale = min(
                1.0, float(self.trainer.current_epoch + 1) / self.hparams["warmup"]
            )
            for pg in optimizer.param_groups:
                pg["lr"] = lr_scale * self.hparams["lr"]

        # update params
        optimizer.step(closure=optimizer_closure)
        optimizer.zero_grad()


class FilterBaseBalanced(FilterBase):
    def __init__(self, hparams):
        super().__init__(hparams)
        """
        Initialise the Lightning Module that can scan over different filter training regimes
        """

    def train_dataloader(self):
        if self.trainset is not None:
            return GeoLoader(self.trainset, batch_size=1, num_workers=1)
        else:
            return None

    def training_step(self, batch, batch_idx):

        # print("Training on device", self.device)
        emb = None if (self.hparams["emb_channels"] == 0) else batch.embedding

        # Handle training towards a subset of the data
        if "subset" in self.hparams["regime"]:
            subset_mask = np.isin(
                batch.edge_index.cpu(), batch.signal_true_edges.unique().cpu()
            ).any(0)
            batch.edge_index = batch.edge_index[:, subset_mask]
            batch.y = batch.y[subset_mask]

        #         print("Starting chunks")
        with torch.no_grad():
            cut_list = []
            for j in range(self.hparams["n_chunks"]):
                # print("Loading chunk", j, "on device", self.device)
                subset_ind = torch.chunk(
                    torch.arange(batch.edge_index.shape[1]), self.hparams["n_chunks"]
                )[j]
                if "ci" in self.hparams["regime"]:
                    output = self(
                        torch.cat(
                            [
                                batch.cell_data[:, : self.hparams["cell_channels"]],
                                batch.x,
                            ],
                            axis=-1,
                        ),
                        batch.edge_index[:, subset_ind],
                        emb,
                    ).squeeze()
                else:
                    output = self(
                        batch.x, batch.edge_index[:, subset_ind], emb
                    ).squeeze()

                cut = torch.sigmoid(output) > self.hparams["filter_cut"]
                cut_list.append(cut)

            cut_list = torch.cat(cut_list)

            num_true, num_false = batch.y.bool().sum(), (~batch.y.bool()).sum()
            true_indices = torch.where(batch.y.bool())[0]
            hard_negatives = cut_list & ~batch.y.bool()
            hard_indices = torch.where(hard_negatives)[0]
            hard_indices = hard_indices[torch.randperm(len(hard_indices))][
                : int(len(true_indices) * self.hparams["ratio"] / 2)
            ]
            easy_indices = torch.where(~batch.y.bool())[0][
                torch.randint(
                    num_false, (int(num_true.item() * self.hparams["ratio"] / 2),)
                )
            ]

            combined_indices = torch.cat([true_indices, hard_indices, easy_indices])

            # Shuffle indices:
            combined_indices = combined_indices[torch.randperm(len(combined_indices))][
                : self.hparams["edges_per_batch"]
            ]
            weight = torch.tensor(self.hparams["weight"])

        if "ci" in self.hparams["regime"]:
            output = self(
                torch.cat(
                    [batch.cell_data[:, : self.hparams["cell_channels"]], batch.x],
                    axis=-1,
                ),
                batch.edge_index[:, combined_indices],
                emb,
            ).squeeze()
        else:
            output = self(batch.x, batch.edge_index[:, combined_indices], emb).squeeze()

        if "weighting" in self.hparams["regime"]:
            manual_weights = batch.weights[combined_indices]
            manual_weights[batch.y[combined_indices] == 0] = 1
        else:
            manual_weights = None

        if "pid" in self.hparams["regime"]:
            y_pid = (
                batch.pid[batch.edge_index[0, combined_indices]]
                == batch.pid[batch.edge_index[1, combined_indices]]
            ) * (batch.pid[batch.edge_index[0, combined_indices]])
            loss = F.binary_cross_entropy_with_logits(
                output, y_pid.float(), weight=manual_weights, pos_weight=weight
            )
        else:
            loss = F.binary_cross_entropy_with_logits(
                output,
                batch.y[combined_indices].float(),
                weight=manual_weights,
                pos_weight=weight,
            )

        self.log("train_loss", loss)
        # print("Returning training loss on device", self.device)
        return loss

    def validation_step(self, batch, batch_idx):

        result = self.shared_evaluation(batch, batch_idx, log=True)
   
        return result

    def test_step(self, batch, batch_idx):

        result = self.shared_evaluation(batch, batch_idx, log=False)

        return result

    def shared_evaluation(self, batch, batch_idx, log=False):

        """
        This method is shared between validation steps and test steps
        """

        # print("Validating on device", self.device)
        
        emb = (
            None if (self.hparams["emb_channels"] == 0) else batch.embedding
        )  # Does this work??

        score_list = []
        val_loss = torch.tensor(0).to(self.device)
        for j in range(self.hparams["n_chunks"]):
            # print("Loading chunk", j, "on device", self.device)
            
            subset_ind = torch.chunk(
                torch.arange(batch.edge_index.shape[1]), self.hparams["n_chunks"]
            )[j]

            if "ci" in self.hparams["regime"]:
                output = self(
                    torch.cat(
                        [batch.cell_data[:, : self.hparams["cell_channels"]], batch.x],
                        axis=-1,
                    ),
                    batch.edge_index[:, subset_ind],
                    emb,
                ).squeeze()
            else:
                output = self(batch.x, batch.edge_index[:, subset_ind], emb).squeeze()

            scores = torch.sigmoid(output)
            score_list.append(scores)

            if "weighting" in self.hparams["regime"]:
                manual_weights = batch.weights[subset_ind]
                manual_weights[batch.y[subset_ind] == 0] = 1
            else:
                manual_weights = None

            if "pid" not in self.hparams["regime"]:
                val_loss = val_loss + F.binary_cross_entropy_with_logits(
                    output, batch.y[subset_ind].float(), weight=manual_weights
                )
            else:
                y_pid = (
                    batch.pid[batch.edge_index[0, subset_ind]]
                    == batch.pid[batch.edge_index[1, subset_ind]]
                )
                val_loss = +F.binary_cross_entropy_with_logits(
                    output, y_pid.float(), weight=manual_weights
                )

        score_list = torch.cat(score_list)
        cut_list = score_list > self.hparams["filter_cut"]

        # Edge filter performance
        edge_positive = cut_list.sum().float()
        if "pid" in self.hparams["regime"]:
            truth = batch.pid[batch.edge_index[0]] == batch.pid[batch.edge_index[1]]
        else:
            truth = batch.y.bool()

        eff, pur, auc = self.get_metrics(truth, score_list)

        if log:
            current_lr = self.optimizers().param_groups[0]["lr"]
            self.log_dict(
                {
                    "eff": eff,
                    "pur": pur,
                    "val_loss": val_loss,
                    "current_lr": current_lr,
                    "auc": auc,
                },
                # sync_dist=True,
                # on_epoch=True
            )
            
        # print("Returning validation loss on device", self.device, )
        return {"loss": val_loss, "preds": score_list, "truth": truth}

class LargeFilterBaseBalanced(FilterBaseBalanced):
    def __init__(self, hparams):
        super().__init__(hparams)
        
    def setup(self, stage):
        # Handle any subset of [train, val, test] data split, assuming that ordering

        input_dirs = [None, None, None]
        input_dirs[: len(self.hparams["datatype_names"])] = [
            os.path.join(self.hparams["input_dir"], datatype)
            for datatype in self.hparams["datatype_names"]
        ]
        
        if "trainset" not in self.__dict__.keys():
            self.trainset, self.valset, self.testset = [
                LargeDataset(input_dir, split, name, self.hparams)
                for split, name, input_dir in zip(self.hparams["datatype_split"], self.hparams["datatype_names"], input_dirs)
            ]
            
        if (
            "logger" in self.__dict__.keys()
            and "_experiment" in self.logger.__dict__.keys()
        ):
            self.logger.experiment.define_metric("val_loss", summary="min")
            self.logger.experiment.define_metric("auc", summary="max")
        
    def train_dataloader(self):
        if self.trainset is not None:
            return GeoLoader(self.trainset, batch_size=1, num_workers=8)
        else:
            return None

    def val_dataloader(self):
        if self.valset is not None:
            num_val_workers = 0 if ("gpus" in self.hparams.keys() and self.hparams["gpus"] > 0) else 8
            return GeoLoader(self.valset, batch_size=1, num_workers=num_val_workers)
        else:
            return None

    def test_dataloader(self):
        if self.testset is not None:
            return GeoLoader(self.testset, batch_size=1, num_workers=0)
        else:
            return None

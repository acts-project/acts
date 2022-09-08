import sys
import os
import copy
import logging
import tracemalloc
import gc
from memory_profiler import profile

from pytorch_lightning.callbacks import Callback
import torch.nn.functional as F
import sklearn.metrics
import matplotlib.pyplot as plt
import torch
import numpy as np
from sklearn.metrics import roc_curve

"""
Class-based Callback inference for integration with Lightning
"""


class GNNTelemetry(Callback):

    """
    This callback contains standardised tests of the performance of a GNN
    """

    def __init__(self):
        super().__init__()
        logging.info("Constructing telemetry callback")

    def on_test_start(self, trainer, pl_module):

        """
        This hook is automatically called when the model is tested after training. The best checkpoint is automatically loaded
        """
        self.preds = []
        self.truth = []

        print("Starting TELEMETRY")

    def on_test_batch_end(
        self, trainer, pl_module, outputs, batch, batch_idx, dataloader_idx
    ):

        """
        Get the relevant outputs from each batch
        """

        self.preds.append(outputs["preds"].cpu())
        self.truth.append(outputs["truth"].cpu())

    def on_test_end(self, trainer, pl_module):

        """
        1. Aggregate all outputs,
        2. Calculate the ROC curve,
        3. Plot ROC curve,
        4. Save plots to PDF 'metrics.pdf'
        """

        metrics = self.calculate_metrics()

        metrics_plots = self.plot_metrics(metrics)

        self.save_metrics(metrics_plots, pl_module.hparams.output_dir)

    def get_eff_pur_metrics(self):

        self.truth = torch.cat(self.truth)
        self.preds = torch.cat(self.preds)

        fpr, eff, score_cuts = roc_curve(self.truth, self.preds)
        pur = 1 - fpr

        eff, pur, score_cuts = (
            eff[score_cuts <= 1],
            pur[score_cuts <= 1],
            score_cuts[score_cuts <= 1],
        )  # Make sure this is nicely plottable!

        return eff, pur, score_cuts

    def calculate_metrics(self):

        eff, pur, score_cuts = self.get_eff_pur_metrics()

        return {
            "eff_plot": {"eff": eff, "score_cuts": score_cuts},
            "pur_plot": {"pur": pur, "score_cuts": score_cuts},
            "auc_plot": {"eff": eff, "pur": pur},
        }

    def make_plot(self, x_val, y_val, x_lab, y_lab, title):

        # Update this to dynamically adapt to number of metrics
        fig, axs = plt.subplots(nrows=1, ncols=1, figsize=(20, 20))
        axs = axs.flatten() if type(axs) is list else [axs]

        axs[0].plot(x_val, y_val)
        axs[0].set_xlabel(x_lab)
        axs[0].set_ylabel(y_lab)
        axs[0].set_title(title)
        plt.tight_layout()

        return fig, axs

    def plot_metrics(self, metrics):

        eff_fig, eff_axs = self.make_plot(
            metrics["eff_plot"]["score_cuts"],
            metrics["eff_plot"]["eff"],
            "cut",
            "Eff",
            "Efficiency vs. cut",
        )
        pur_fig, pur_axs = self.make_plot(
            metrics["pur_plot"]["score_cuts"],
            metrics["pur_plot"]["pur"],
            "cut",
            "Pur",
            "Purity vs. cut",
        )

        auc_fig, auc_axs = self.make_plot(
            metrics["auc_plot"]["eff"],
            metrics["auc_plot"]["pur"],
            "Eff",
            "Pur",
            "Purity vs. Eff",
        )

        return {
            "eff_plot": [eff_fig, eff_axs],
            "pur_plot": [pur_fig, pur_axs],
            "auc_plot": [auc_fig, auc_axs],
        }

    def save_metrics(self, metrics_plots, output_dir):

        os.makedirs(output_dir, exist_ok=True)

        for metric, (fig, axs) in metrics_plots.items():
            fig.savefig(os.path.join(output_dir, f"metrics_{metric}.pdf"), format="pdf")


class GNNBuilder(Callback):
    """Callback handling filter inference for later stages.

    This callback is used to apply a trained filter model to the dataset of a LightningModule.
    The data structure is preloaded in the model, as training, validation and testing sets.
    Intended usage: run training and examine the telemetry to decide on the hyperparameters (e.g. r_test) that
    lead to desired efficiency-purity tradeoff. Then set these hyperparameters in the pipeline configuration and run
    with the --inference flag. Otherwise, to just run straight through automatically, train with this callback included.

    """

    def __init__(self):
        self.output_dir = None
        self.overwrite = False

    def on_test_end(self, trainer, pl_module):

        print("Testing finished, running inference to build graphs...")

        datasets = self.prepare_datastructure(pl_module)

        total_length = sum([len(dataset) for dataset in datasets.values()])

        pl_module.eval()
        with torch.no_grad():
            batch_incr = 0
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
                        self.construct_downstream(batch_to_save, pl_module, datatype)

                    batch_incr += 1

    def prepare_datastructure(self, pl_module):
        # Prep the directory to produce inference data to
        self.output_dir = pl_module.hparams.output_dir
        self.datatypes = ["train", "val", "test"]

        os.makedirs(self.output_dir, exist_ok=True)
        [
            os.makedirs(os.path.join(self.output_dir, datatype), exist_ok=True)
            for datatype in self.datatypes
        ]

        # Set overwrite setting if it is in config
        self.overwrite = (
            pl_module.hparams.overwrite if "overwrite" in pl_module.hparams else False
        )

        # By default, the set of examples propagated through the pipeline will be train+val+test set
        datasets = {
            "train": pl_module.trainset,
            "val": pl_module.valset,
            "test": pl_module.testset,
        }

        return datasets

    def construct_downstream(self, batch, pl_module, datatype):

        output = pl_module(
            pl_module.get_input_data(batch),
            torch.cat([batch.edge_index, batch.edge_index.flip(0)], dim=-1),
        ).squeeze()
        batch.scores = torch.sigmoid(output)

        self.save_downstream(batch, pl_module, datatype)

    def save_downstream(self, batch, pl_module, datatype):

        with open(
            os.path.join(self.output_dir, datatype, os.path.basename(batch.event_file)), "wb"
        ) as pickle_file:
            torch.save(batch, pickle_file)

        logging.info("Saved event {}".format(batch.event_file[-4:]))

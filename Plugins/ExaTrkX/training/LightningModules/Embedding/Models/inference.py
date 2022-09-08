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

from ..utils import build_edges, graph_intersection

"""
Class-based Callback inference for integration with Lightning
"""


class EmbeddingTelemetry(Callback):

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
        self.truth_graph = []
        self.pt_true_pos = []
        self.pt_true = []
        self.distances = []

        self.hparams = pl_module.hparams

    def on_test_batch_end(
        self, trainer, pl_module, outputs, batch, batch_idx, dataloader_idx
    ):

        """
        Get the relevant outputs from each batch
        """

        true_positives = outputs["preds"][:, outputs["truth"]]
        true = outputs["truth_graph"]

        if "pt" in batch:
            pts = batch.pt
            self.pt_true_pos.append(pts[true_positives].cpu())
            self.pt_true.append(pts[true].cpu())

        self.truth.append(outputs["truth"].cpu())
        self.distances.append(outputs["distances"].cpu())
        self.truth_graph.append(outputs["truth_graph"].cpu())

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

    def get_pt_metrics(self):

        pt_true_pos = np.concatenate(self.pt_true_pos, axis=1)
        pt_true = np.concatenate(self.pt_true, axis=1)

        pt_true_pos_av = (pt_true_pos[0] + pt_true_pos[1]) / 2
        pt_true_av = (pt_true[0] + pt_true[1]) / 2

        #         bins = np.arange(pl_module.hparams["pt_min"], np.ceil(pt_true_av.max()), 0.5)
        #         bins = np.logspace(np.log(np.floor(pt_true_av.min())), np.log(np.ceil(pt_true_av.max())), 10)
        bins = np.logspace(0, 1.5, 10)
        centers = [(bins[i] + bins[i + 1]) / 2 for i in range(len(bins) - 1)]

        tp_hist = np.histogram(pt_true_pos_av, bins=bins)[0]
        t_hist = np.histogram(pt_true_av, bins=bins)[0]
        ratio_hist = tp_hist / t_hist

        return centers, ratio_hist

    def get_eff_pur_metrics(self):

        self.distances = torch.cat(self.distances)
        self.truth = torch.cat(self.truth)
        self.truth_graph = torch.cat(self.truth_graph, axis=1)

        r_cuts = np.arange(0.01, self.hparams["r_test"], 0.01)

        #print(r_cuts)

        positives = np.array([(self.distances < r_cut**2).sum() for r_cut in r_cuts])
        true_positives = np.array(
            [self.truth[self.distances < r_cut**2].sum() for r_cut in r_cuts]
        )

        #print(positives, true_positives)
        eff = true_positives / self.truth_graph.shape[1]
        pur = true_positives / positives

        return eff, pur, r_cuts

    def calculate_metrics(self):

        centers, ratio_hist = self.get_pt_metrics()

        eff, pur, r_cuts = self.get_eff_pur_metrics()

        return {
            "pt_plot": {"centers": centers, "ratio_hist": ratio_hist},
            "eff_plot": {"eff": eff, "r_cuts": r_cuts},
            "pur_plot": {"pur": pur, "r_cuts": r_cuts},
        }

    def make_plot(self, x_val, y_val, x_lab, y_lab, title):

        # Update this to dynamically adapt to number of metrics
        fig, axs = plt.subplots(nrows=1, ncols=1, figsize=(10, 10))
        axs = axs.flatten() if type(axs) is list else [axs]

        axs[0].plot(x_val, y_val)
        axs[0].set_xlabel(x_lab)
        axs[0].set_ylabel(y_lab)
        axs[0].set_title(title)
        plt.tight_layout()

        return fig, axs

    def plot_metrics(self, metrics):

        centers, ratio_hist = (
            metrics["pt_plot"]["centers"],
            metrics["pt_plot"]["ratio_hist"],
        )
        pt_fig, pt_axs = self.make_plot(
            centers, ratio_hist, "pT (GeV)", "Efficiency", "Metric Learning Efficiency"
        )

        eff_fig, eff_axs = self.make_plot(
            metrics["eff_plot"]["r_cuts"],
            metrics["eff_plot"]["eff"],
            "radius",
            "Eff",
            "Efficiency vs. radius",
        )
        pur_fig, pur_axs = self.make_plot(
            metrics["pur_plot"]["r_cuts"],
            metrics["pur_plot"]["pur"],
            "radius",
            "Pur",
            "Purity vs. radius",
        )
        roc_fig, roc_axs = self.make_plot(
            metrics["pur_plot"]["pur"],
            metrics["eff_plot"]["eff"],
            "Pur",
            "Eff",
            "ROC Curve",
        )

        return {
            "pt_plot": [pt_fig, pt_axs],
            "eff_plot": [eff_fig, eff_axs],
            "pur_plot": [pur_fig, pur_axs],
            "roc_plot": [roc_fig, roc_axs],
        }

    def save_metrics(self, metrics_plots, output_dir):

        os.makedirs(output_dir, exist_ok=True)

        for metric, (fig, axs) in metrics_plots.items():
            fig.savefig(os.path.join(output_dir, f"metrics_{metric}.pdf"), format="pdf")


class EmbeddingBuilder(Callback):
    """Callback handling embedding inference for later stages.

    This callback is used to apply a trained embedding model to the dataset of a LightningModule.
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

    def get_truth(self, batch, e_spatial, e_bidir):

        e_spatial_easy_fake = e_spatial[
            :, batch.pid[e_spatial[0]] != batch.pid[e_spatial[1]]
        ]
        y_cluster_easy_fake = torch.zeros(e_spatial_easy_fake.shape[1])

        e_spatial_ambiguous = e_spatial[
            :, batch.pid[e_spatial[0]] == batch.pid[e_spatial[1]]
        ]
        e_spatial_ambiguous, y_cluster_ambiguous = graph_intersection(
            e_spatial_ambiguous, e_bidir
        )

        e_spatial = torch.cat([e_spatial_easy_fake.cpu(), e_spatial_ambiguous], dim=-1)
        y_cluster = torch.cat([y_cluster_easy_fake, y_cluster_ambiguous])

        return e_spatial, y_cluster

    def construct_downstream(self, batch, pl_module, datatype):

        input_data = pl_module.get_input_data(batch)

        spatial = pl_module(input_data)

        # Make truth bidirectional
        e_bidir = torch.cat(
            [batch.signal_true_edges, batch.signal_true_edges.flip(0)],
            axis=-1,
        )

        # Build the radius graph with radius < r_test
        e_spatial = build_edges(
            spatial, spatial, indices=None, r_max=pl_module.hparams.r_test, k_max=1000
        )  # This step should remove reliance on r_val, and instead compute an r_build based on the EXACT r required to reach target eff/pur

        # Arbitrary ordering to remove half of the duplicate edges
        R_dist = torch.sqrt(batch.x[:, 0] ** 2 + batch.x[:, 2] ** 2)
        e_spatial = e_spatial[:, (R_dist[e_spatial[0]] <= R_dist[e_spatial[1]])]

        e_spatial, y_cluster = self.get_truth(batch, e_spatial, e_bidir)

        # Re-introduce random direction, to avoid training bias
        random_flip = torch.randint(2, (e_spatial.shape[1],)).bool()
        e_spatial[0, random_flip], e_spatial[1, random_flip] = (
            e_spatial[1, random_flip],
            e_spatial[0, random_flip],
        )

        batch.edge_index = e_spatial
        batch.y = y_cluster
        batch.y_pid = batch.pid[batch.edge_index[0]] == batch.pid[batch.edge_index[1]]

        self.save_downstream(batch, pl_module, datatype)

    def save_downstream(self, batch, pl_module, datatype):

        with open(
            os.path.join(self.output_dir, datatype, os.path.basename(batch.event_file)), "wb"
        ) as pickle_file:
            torch.save(batch, pickle_file)

        logging.info("Saved event {}".format(batch.event_file[-4:]))

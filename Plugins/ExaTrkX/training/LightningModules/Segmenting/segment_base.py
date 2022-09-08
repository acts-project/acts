"""A base class for the segmentation data module.

The segmentation data module can be used as a stand-alone processing module for converted classified (i.e. scored) graphs into track candidates. It fills a gap between the initial GNN (which has some tracking telemetry, but only scores graphs), and any later stages that can train on segments (but take graph scores and form segments on-the-fly).

Example:
    This model can be run with `traintrack` as with other stages. One can then experiment with the segments that are produced.

"""

# 3rd party imports
from pytorch_lightning import LightningDataModule


class SegmentBase(LightningDataModule):
    def __init__(self, hparams):
        super().__init__()
        self.save_hyperparameters(hparams)

        self.input_dir = self.hparams["input_dir"]
        self.output_dir = self.hparams["output_dir"]
        self.n_files = self.hparams["n_files"]
        self.edge_cut = (
            0.5 if "edge_cut" not in self.hparams else self.hparams["edge_cut"]
        )

        self.n_tasks = (
            1 if "n_tasks" not in self.hparams else self.hparams["n_tasks"]
        )
        self.task = 0 if "task" not in self.hparams else self.hparams["task"]
        self.n_workers = (
            self.hparams["n_workers"]
            if "n_workers" in self.hparams
            else len(os.sched_getaffinity(0))
        )

        self.show_progress = (
            self.hparams["show_progress"] if "show_progress" in self.hparams else True
        )

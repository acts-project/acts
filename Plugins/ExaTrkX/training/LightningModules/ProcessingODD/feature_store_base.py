# 3rd party imports
from pytorch_lightning import LightningDataModule


class FeatureStoreBase(LightningDataModule):
    def __init__(self, hparams):
        super().__init__()
        self.save_hyperparameters(hparams)

        self.input_dir = self.hparams["input_dir"]
        self.output_dir = self.hparams["output_dir"]
        self.n_files = self.hparams["n_files"]

        self.n_tasks = self.hparams["n_tasks"]
        self.task = 0 if "task" not in self.hparams else self.hparams["task"]
        self.n_workers = (
            self.hparams["n_workers"]
            if "n_workers" in self.hparams
            else len(os.sched_getaffinity(0))
        )
        self.build_weights = (
            self.hparams["build_weights"] if "build_weights" in self.hparams else True
        )
        self.show_progress = (
            self.hparams["show_progress"] if "show_progress" in self.hparams else True
        )

from tqdm import tqdm
from tqdm.contrib.concurrent import process_map

from ..segment_base import SegmentBase
from ..utils.segmentation_utils import label_graph

import os
import numpy as np

from functools import partial


class TrackMLSegment(SegmentBase):

    """
    The segmentation data module specific to the TrackML pipeline
    """

    def __init__(self, hparams):
        super().__init__(hparams)

        """Init method for class.
        
        Args:
            params (type): Description of params.
            
        """

    def prepare_data(self):

        all_files = []
        
        for subdir in ["train", "test", "val"]:
            all_files += [
                os.path.join(self.hparams["input_dir"], subdir, file)
                for file in os.listdir(os.path.join(self.hparams["input_dir"], subdir))
            ]
            
        all_files = all_files[:self.n_files]
            
        all_files = np.array_split(all_files, self.n_tasks)[self.task]

        os.makedirs(self.output_dir, exist_ok=True)
        print("Writing outputs to " + self.output_dir)

        process_func = partial(label_graph, **self.hparams)
        process_map(process_func, all_files, max_workers=self.n_workers)

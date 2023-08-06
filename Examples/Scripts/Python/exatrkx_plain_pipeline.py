#!/usr/bin/env python3

from multiprocessing import Process, Queue
from pathlib import Path
import argparse
import os

import numpy as np


def load_pytorch_graph(queue, filename):
    import torch

    graph = torch.load(filename)

    graphdict = {}
    for k in graph.keys:
        try:
            graphdict[k] = graph[k].detach().cpu().numpy()
        except:
            continue

    queue.put(graphdict)


def run_pipeline(data, torchscript_dir):
    import acts
    import acts.examples
    
    # Make input data
    scaling = ScalingDict(
        {
            "r": 1000.0,
            "phi": 3.14159,
            "z": 3000.0,
        }
    )

    input_tensor = np.vstack(
        [
            data[k] / scaling[k]
            for k in ["r", "phi", "z", "cell_count", "cell_val", "lx", "ly"]
        ]
    ).T

    # Make metric hook
    hook = acts.examples.TorchTruthGraphMetricsHook(
        data["track_edges"].T.flatten(), acts.logging.DEBUG
    )

    # Make stages
    emb = acts.examples.TorchMetricLearning(
        acts.logging.VERBOSE,
        embeddingDim=8,
        knnVal=100,
        numFeatures=7,
        rVal=0.2,
        modelPath=torchscript_dir / "embedding.pt",
    )

    flt = acts.examples.TorchEdgeClassifier(
        acts.logging.VERBOSE,
        numFeatures=3,
        undirected=False,
        modelPath=torchscript_dir / "filter.pt",
    )

    gnn = acts.examples.TorchEdgeClassifier(
        acts.logging.VERBOSE,
        undirected=True,
        numFeatures=3,
        modelPath=torchscript_dir / "gnn.pt",
    )

    trk = acts.examples.BoostTrackBuilding(acts.logging.VERBOSE)
    
    # Run pipeline
    pipeline = acts.examples.Pipeline(emb, [flt,gnn], trk, acts.logging.VERBOSE)

    return pipeline.run(
        input_tensor.flatten().tolist(), np.arange(input_tensor.shape[0]).tolist(), hook
    )


class ScalingDict(dict):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def __missing__(self, key):
        return 1.0


def main():

    parser = argparse.ArgumentParser()
    parser.add_argument("torchscript_dir", type=str)
    parser.add_argument("pyg_file", type=str)
    args = vars(parser.parse_args())

    assert os.path.exists(args["torchscript_dir"])
    assert os.path.exists(args["pyg_file"])

    # Load *.pyg in different process to avoid conflicting torch versions for python and C++
    q = Queue()
    p = Process(target=load_pytorch_graph, args=(q, args["pyg_file"]))
    p.start()
    graph = q.get()
    p.join()

    run_pipeline(graph, Path(args["torchscript_dir"]))



if __name__ == "__main__":
    main()

"""
TODO:

- Change manual scipy sparse conversion to PyG version for brevity
"""

import os
import logging

import torch
import scipy.sparse.csgraph as scigraph
import scipy.sparse as sp
import numpy as np
from torch_geometric.utils import to_scipy_sparse_matrix


def label_graph(
    input_file: str, output_dir: str, edge_cut: float = 0.5, **kwargs
) -> None:

    """Loads an input_file and outputs a segmented (i.e. labelled) graph.

    Args:
        input_file: Location of the input graph (a torch pickled file containing a Pytorch Geometric data object).
        edge_cut: The minimum score for an edge to become part of a segment

    """

    #try:
    if True:
        output_file = os.path.join(output_dir, os.path.split(input_file)[-1])

        overwrite=True
        if not os.path.exists(output_file) or overwrite:

            logging.info("Preparing event {}".format(output_file))
            graph = torch.load(input_file, map_location="cpu")

            # apply cut
            passing_edges = graph.edge_index[:, graph.scores[:graph.edge_index.shape[1]] > edge_cut]

            # get connected components
            sparse_edges = sp.coo_matrix(
                (np.ones(passing_edges.shape[1]), passing_edges.cpu().numpy()),
                shape=(len(graph.x), len(graph.x)),
            )
            connected_components = scigraph.connected_components(sparse_edges)[1]

            # attach labels to data
            graph.labels = connected_components

            with open(output_file, "wb") as pickle_file:
                torch.save(graph, pickle_file)

        else:
            logging.info("{} already exists".format(output_file))

    #except Exception as inst:
        #print("File:", input_file, "had exception", inst)

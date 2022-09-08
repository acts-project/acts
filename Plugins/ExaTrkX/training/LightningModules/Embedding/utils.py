import os
import logging

import torch
from torch.utils.data import random_split
from torch import nn
import scipy as sp
import numpy as np

from alive_progress import alive_bar


"""
Ideally, we would be using FRNN and the GPU. But in the case of a user not having a GPU, or not having FRNN, we import FAISS as the 
nearest neighbor library
"""

import faiss
import faiss.contrib.torch_utils

try:
    import frnn

    FRNN_AVAILABLE = True
except ImportError:
    FRNN_AVAILABLE = False

if torch.cuda.is_available():
    device = "cuda"
else:
    device = "cpu"
    FRNN_AVAILABLE = False

#FRNN_AVAILABLE = False

def load_dataset(
    input_dir,
    num,
    pt_background_cut,
    pt_signal_cut,
    nhits,
    primary_only,
    true_edges,
    noise,
    **kwargs
):
    if input_dir is not None:
        all_events = os.listdir(input_dir)
        all_events = sorted([os.path.join(input_dir, event) for event in all_events])
        loaded_events = []

        disable_bar = not ("ALIVE_BAR" in os.environ and bool(os.environ["ALIVE_BAR"]))
        with alive_bar(num, disable=disable_bar, title="Loading") as bar:
            for event in all_events[:num]:
                try:
                    loaded_event = torch.load(event, map_location=torch.device("cpu"))
                    loaded_event.event_file = event
                    loaded_events.append(loaded_event)
                except:
                    logging.info("Corrupted event file: {}".format(event))

                bar()

        loaded_events = select_data(
            loaded_events,
            pt_background_cut,
            pt_signal_cut,
            nhits,
            primary_only,
            true_edges,
            noise,
        )
        return loaded_events
    else:
        return None


def split_datasets(
    input_dir="",
    train_split=[100, 10, 10],
    pt_background_cut=0,
    pt_signal_cut=0,
    nhits=0,
    primary_only=False,
    true_edges=None,
    noise=True,
    seed=1,
    **kwargs
):
    """
    Prepare the random Train, Val, Test split, using a seed for reproducibility. Seed should be
    changed across final varied runs, but can be left as default for experimentation.
    """

    torch.manual_seed(seed)
    loaded_events = load_dataset(
        input_dir,
        sum(train_split),
        pt_background_cut,
        pt_signal_cut,
        nhits,
        primary_only,
        true_edges,
        noise,
    )
    train_events, val_events, test_events = random_split(loaded_events, train_split)

    return train_events, val_events, test_events


def get_edge_subset(edges, mask_where, inverse_mask):

    included_edges_mask = np.isin(edges, mask_where).all(0)
    included_edges = edges[:, included_edges_mask]
    included_edges = inverse_mask[included_edges]

    return included_edges, included_edges_mask


def select_data(
    events, pt_background_cut, pt_signal_cut, nhits_min, primary_only, true_edges, noise
):
    # Handle event in batched form
    if type(events) is not list:
        events = [events]

    # NOTE: Cutting background by pT BY DEFINITION removes noise
    if pt_background_cut > 0 or not noise:
        disable_bar = not ("ALIVE_BAR" in os.environ and bool(os.environ["ALIVE_BAR"]))
        with alive_bar(len(events), disable=disable_bar, title="Processing I") as bar:
            for event in events:

                pt_mask = (event.pt > pt_background_cut) & (event.pid == event.pid) & (event.pid != 0)
                pt_where = torch.where(pt_mask)[0]

                inverse_mask = torch.zeros(pt_where.max() + 1).long()
                inverse_mask[pt_where] = torch.arange(len(pt_where))

                event[true_edges], edge_mask = get_edge_subset(
                    event[true_edges], pt_where, inverse_mask
                )

                node_features = ["cell_data", "x", "hid", "pid", "pt", "nhits", "primary"]
                for feature in node_features:
                    if feature in event.keys:
                        event[feature] = event[feature][pt_mask]

                bar()

    disable_bar = not ("ALIVE_BAR" in os.environ and bool(os.environ["ALIVE_BAR"]))
    with alive_bar(len(events), disable=disable_bar, title="Processing II") as bar:
        for event in events:
            event.signal_true_edges = event[true_edges]
            edge_subset = torch.ones(event.signal_true_edges.shape[1]).bool()

            if "pt" in event.keys:
                edge_subset &= (event.pt[event[true_edges]] > pt_signal_cut).all(0)

            if "primary" in event.keys:
                edge_subset &= (event.nhits[event[true_edges]] >= nhits_min).all(0)

            if "nhits" in event.keys:
                edge_subset &= ((event.primary[event[true_edges]].bool().all(0) | (not primary_only)))

            event.signal_true_edges = event.signal_true_edges[:, edge_subset]

            bar()

    return events


def reset_edge_id(subset, graph):
    subset_ind = np.where(subset)[0]
    filler = -np.ones((graph.max() + 1,))
    filler[subset_ind] = np.arange(len(subset_ind))
    graph = torch.from_numpy(filler[graph]).long()
    exist_edges = (graph[0] >= 0) & (graph[1] >= 0)
    graph = graph[:, exist_edges]

    return graph, exist_edges


def graph_intersection(
    pred_graph, truth_graph, using_weights=False, weights_bidir=None
):

    array_size = max(pred_graph.max().item(), truth_graph.max().item()) + 1

    if torch.is_tensor(pred_graph):
        l1 = pred_graph.cpu().numpy()
    else:
        l1 = pred_graph
    if torch.is_tensor(truth_graph):
        l2 = truth_graph.cpu().numpy()
    else:
        l2 = truth_graph
    e_1 = sp.sparse.coo_matrix(
        (np.ones(l1.shape[1]), l1), shape=(array_size, array_size)
    ).tocsr()
    e_2 = sp.sparse.coo_matrix(
        (np.ones(l2.shape[1]), l2), shape=(array_size, array_size)
    ).tocsr()
    del l1

    e_intersection = e_1.multiply(e_2) - ((e_1 - e_2) > 0)
    del e_1
    del e_2

    if using_weights:
        weights_list = weights_bidir.cpu().numpy()
        weights_sparse = sp.sparse.coo_matrix(
            (weights_list, l2), shape=(array_size, array_size)
        ).tocsr()
        del weights_list
        del l2
        new_weights = weights_sparse[e_intersection.astype("bool")]
        del weights_sparse
        new_weights = torch.from_numpy(np.array(new_weights)[0])

    e_intersection = e_intersection.tocoo()
    new_pred_graph = torch.from_numpy(
        np.vstack([e_intersection.row, e_intersection.col])
    ).long()  # .to(device)
    y = torch.from_numpy(e_intersection.data > 0)  # .to(device)
    del e_intersection

    if using_weights:
        return new_pred_graph, y, new_weights
    else:
        return new_pred_graph, y


def build_edges(
    query, database, indices=None, r_max=1.0, k_max=10, return_indices=False
):
    """
    NOTE: These KNN/FRNN algorithms return the distances**2. Therefore we need to be careful when comparing them to the target distances (r_val, r_test), and to the margin parameter (which is L1 distance)
    """

    if FRNN_AVAILABLE:

        Dsq, I, nn, grid = frnn.frnn_grid_points(
            points1=query.unsqueeze(0),
            points2=database.unsqueeze(0),
            lengths1=None,
            lengths2=None,
            K=k_max,
            r=r_max,
            grid=None,
            return_nn=False,
            return_sorted=True,
        )

        I = I.squeeze().int()
        ind = torch.Tensor.repeat(
            torch.arange(I.shape[0], device=device), (I.shape[1], 1), 1
        ).T.int()
        positive_idxs = I >= 0
        edge_list = torch.stack([ind[positive_idxs], I[positive_idxs]]).long()

    else:

        if device == "cuda":
            res = faiss.StandardGpuResources()
            Dsq, I = faiss.knn_gpu(res=res, xq=query, xb=database, k=k_max)
        elif device == "cpu":
            index = faiss.IndexFlatL2(database.shape[1])
            index.add(database)
            Dsq, I = index.search(query, k_max)

        ind = torch.Tensor.repeat(
            torch.arange(I.shape[0], device=device), (I.shape[1], 1), 1
        ).T.int()

        edge_list = torch.stack([ind[Dsq <= r_max**2], I[Dsq <= r_max**2]])

    # Reset indices subset to correct global index
    if indices is not None:
        edge_list[0] = indices[edge_list[0]]

    # Remove self-loops
    edge_list = edge_list[:, edge_list[0] != edge_list[1]]

    if return_indices:
        return edge_list, Dsq, I, ind
    else:
        return edge_list


def build_knn(spatial, k):

    if device == "cuda":
        res = faiss.StandardGpuResources()
        _, I = faiss.knn_gpu(res=res, xq=spatial, xb=spatial, k=k_max)
    elif device == "cpu":
        index = faiss.IndexFlatL2(spatial.shape[1])
        index.add(spatial)
        _, I = index.search(spatial, k_max)

    ind = torch.Tensor.repeat(
        torch.arange(I.shape[0], device=device), (I.shape[1], 1), 1
    ).T
    edge_list = torch.stack([ind, I])

    # Remove self-loops
    edge_list = edge_list[:, edge_list[0] != edge_list[1]]

    return edge_list


def get_best_run(run_label, wandb_save_dir):
    for (root_dir, dirs, files) in os.walk(wandb_save_dir + "/wandb"):
        if run_label in dirs:
            run_root = root_dir

    best_run_base = os.path.join(run_root, run_label, "checkpoints")
    best_run = os.listdir(best_run_base)
    best_run_path = os.path.join(best_run_base, best_run[0])

    return best_run_path


# -------------------------- Performance Evaluation -------------------


def embedding_model_evaluation(model, trainer, fom="eff", fixed_value=0.96):

    # Seed solver with one batch, then run on full test dataset
    sol = root(
        evaluate_set_root,
        args=(model, trainer, fixed_value, fom),
        x0=0.9,
        x1=1.2,
        xtol=0.001,
    )
    print("Seed solver complete, radius:", sol.root)

    # Return ( (efficiency, purity), radius_size)
    return evaluate_set_metrics(sol.root, model, trainer), sol.root


def evaluate_set_root(r, model, trainer, goal=0.96, fom="eff"):
    eff, pur = evaluate_set_metrics(r, model, trainer)

    if fom == "eff":
        return eff - goal

    elif fom == "pur":
        return pur - goal


def get_metrics(test_results, model):

    ps = [len(result["truth"]) for result in test_results]
    ts = [result["truth_graph"].shape[1] for result in test_results]
    tps = [result["truth"].sum() for result in test_results]

    efficiencies = [tp / t for (t, tp) in zip(ts, tps)]
    purities = [tp / p for (p, tp) in zip(ps, tps)]

    mean_efficiency = np.mean(efficiencies)
    mean_purity = np.mean(purities)

    return mean_efficiency, mean_purity


def evaluate_set_metrics(r_test, model, trainer):

    model.hparams.r_test = r_test
    test_results = trainer.test(ckpt_path=None)

    mean_efficiency, mean_purity = get_metrics(test_results, model)

    print(mean_purity, mean_efficiency)

    return mean_efficiency, mean_purity


# ------------------------- Convenience Utilities ---------------------------


def make_mlp(
    input_size,
    sizes,
    hidden_activation="ReLU",
    output_activation="ReLU",
    layer_norm=False,
):
    """Construct an MLP with specified fully-connected layers."""
    hidden_activation = getattr(nn, hidden_activation)
    if output_activation is not None:
        output_activation = getattr(nn, output_activation)
    layers = []
    n_layers = len(sizes)
    sizes = [input_size] + sizes
    # Hidden layers
    for i in range(n_layers - 1):
        layers.append(nn.Linear(sizes[i], sizes[i + 1]))
        if layer_norm:
            layers.append(nn.LayerNorm(sizes[i + 1]))
        layers.append(hidden_activation())
    # Final layer
    layers.append(nn.Linear(sizes[-2], sizes[-1]))
    if output_activation is not None:
        if layer_norm:
            layers.append(nn.LayerNorm(sizes[-1]))
        layers.append(output_activation())
    return nn.Sequential(*layers)

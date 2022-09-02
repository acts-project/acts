# Exa.TrkX plugin

This plugin contains a track finding module based on Graph Neural Networks (GNNs) which is developed by the [Exa.TrkX](https://exatrkx.github.io/) team. Build instructions and dependencies can be found in the [README](https://github.com/acts-project/acts/blob/main/Plugins/ExaTrkX/README.md) of the plugin.

## Model

The Exa.TrkX is a three-stage GNN-based algorithm:

1) **Graph building** (*Embedding stage*): This is done via metric learning approach. A neural network tries to learn an mapping that minimizes the distance between points of the same track in the embedded space. In this embedded space is then a graph built using a fixed nearest-neighbor search.
2) **Graph size reduction** (*Filter stage*): The resulting graph is typically to for processing with a graph neural network. Therefore a simple classifier reduces the graph size further.
3) **Edge classification** (*GNN stage*): Finally, a GNN is used to find the edges in the Graph that belong to tracks.

Finally, the track candidates must be extracted by a algorithm that finds the `weakly connected components`. For the implementation, see below.

## Implementation

### Neural network backend

Currently there are two backends available:

- The **TorchScript** backend requires besides `libtorch` also `TorchScatter` that implements the scatter operators used in the graph neural network. 
- The **ONNX** backend is currently not as well maintained as the TorchScript backend. The main reason is, that there are scatterin operators missing in ONNX, so that the graph neural network cannot be exported correctely in the `.onnx` format. The implementation also uses `cugraph` for graph building instead of `Boost.Graph`.

### Graph building

Both backends use currently different libraries for graph building.

- The TorchScript backend uses `Boost.Graph` for graph building.
- The ONNX backend uses the `cugraph` library for graph building

## Ressources

* Talk by *Daniel Murnane* at the [Connecting the Dots 2020](https://indico.cern.ch/event/831165/contributions/3717124/attachments/2024241/3385587/GNNs_for_Track_Finding.pdf)
* Talk by *Daniel Murnane* at the [vCHEP 2021](https://indico.cern.ch/event/948465/contributions/4323753/attachments/2246789/3810686/Physics%20and%20Computing%20Performance%20of%20the%20ExaTrkX%20TrackML%20Pipeline.pdf)
* Talk by *Alina Lazar* at the [ACAT 2021](https://indico.cern.ch/event/855454/contributions/4605079/attachments/2357191/4022841/ExaTrkX%20Inference%20-%20ACAT21%20v7.pdf)
* Talk by *Benjamin Huth* at the [ICHEP 2022](https://agenda.infn.it/event/28874/contributions/169199/attachments/94163/128944/slides_benjamin_huth_exatkrkx_acts.pdf)

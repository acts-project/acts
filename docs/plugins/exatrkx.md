(exatrkxplugin)=
# Exa.TrkX plugin

This plugin contains a track finding module based on Graph Neural Networks (GNNs) which is developed by the [Exa.TrkX](https://exatrkx.github.io/) team. Build instructions and dependencies can be found in the [README](https://github.com/acts-project/acts/blob/main/Plugins/ExaTrkX/README.md) of the plugin.

## Model

The Exa.TrkX pipeline is a three-stage GNN-based algorithm:

1) **Graph building** (*Embedding stage*): This is done via metric learning approach. A neural network tries to learn a mapping that minimizes the distance between points of the same track in the embedded space. In this embedded space then a graph is built using a fixed nearest-neighbor search.
2) **Graph size reduction** (*Filter stage*): The resulting graph is typically to large for processing with a GNN. Therefore a binary classifier reduces the graph size further and removes easy to identify not matching edges.
3) **Edge classification** (*GNN stage*): Finally, a GNN is used to find the edges in the Graph that belong to tracks. This is done by scoring each edge with a number between 0 and 1, and a corresponding *edge cut*.

Finally, the track candidates must be extracted by a algorithm that finds the *weakly connected components* in the graph depending on the *edge cut*.

## Implementation

### Neural network backend

Currently there are two backends available:

- The **TorchScript** backend requires besides *libtorch* also *TorchScatter* that implements the scatter operators used in the graph neural network. The corresponding class is {class}`Acts::ExaTrkXTrackFindingTorch`.
- The **ONNX** backend is currently not as well maintained as the TorchScript backend. The main reason is, that some scattering operators are not available in ONNX, so that the graph neural network cannot be exported correctely in the `.onnx` format. The corresponding class is {class}`Acts::ExaTrkXTrackFindingOnnx`.

### Graph building

Both backends use currently different libraries for graph building.

- The TorchScript backend uses *Boost.Graph* for graph building.
- The ONNX backend uses the *cugraph* library for graph building

## API and examples integration

The interface of the backends is defined by {class}`Acts::ExaTrkXTrackFindingBase` which is also used by `ActsExamples::TrackFindingAlgorithmExaTrkX`. The inference can be called with the `getTracks` function:

```{doxygenclass} Acts::ExaTrkXTrackFindingBase
---
outline:
members: getTracks
---
```

This function takes the the input data as a `std::vector<double>` (e.g., a flattened $N \times 3$ array in cylindric coordinates like $[r_0, \varphi_0, z_0, r_1, \dots \varphi_N, z_N]$), as well as some corresponding spacepoint ids as a `std::vector<int>`. It then fills a `std::vector<std::vector<int>>` with the found track candidates using the provided spacepoint ids. Logging and timing measurements can be enabled with the remaining arguments. The hyperparameters of the models are defined in the `Config` member structs.

:::{note}
Any kind of preprocessing (scaling, ...) of the input values must be done before passing them to the inference.
:::

See [here](https://github.com/acts-project/acts/blob/main/Examples/Scripts/Python/exatrkx.py) for the corresponding python example.

## Ressources

* Talk by *Daniel Murnane* at the [Connecting the Dots 2020](https://indico.cern.ch/event/831165/contributions/3717124/attachments/2024241/3385587/GNNs_for_Track_Finding.pdf)
* Talk by *Daniel Murnane* at the [vCHEP 2021](https://indico.cern.ch/event/948465/contributions/4323753/attachments/2246789/3810686/Physics%20and%20Computing%20Performance%20of%20the%20ExaTrkX%20TrackML%20Pipeline.pdf)
* Talk by *Alina Lazar* at the [ACAT 2021](https://indico.cern.ch/event/855454/contributions/4605079/attachments/2357191/4022841/ExaTrkX%20Inference%20-%20ACAT21%20v7.pdf)
* Talk by *Benjamin Huth* at the [ICHEP 2022](https://agenda.infn.it/event/28874/contributions/169199/attachments/94163/128944/slides_benjamin_huth_exatkrkx_acts.pdf)

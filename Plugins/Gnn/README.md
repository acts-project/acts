# GNN Plugin

This plugin contains a track finding module based on Graph Neural Networks (GNNs). More background information can be found in the [documentation](https://acts.readthedocs.io/en/latest/plugins/gnn.html).

## Building

To build the plugin, enable the appropriate CMake options:

```bash
cmake -B <build> -S <source> \
  -D ACTS_BUILD_PLUGIN_GNN=ON \
  -D ACTS_GNN_ENABLE_TORCH=ON/OFF \
  -D ACTS_GNN_ENABLE_ONNX=ON/OFF \
  -D ACTS_BUILD_EXAMPLES_GNN=ON \
  -D ACTS_BUILD_EXAMPLES_PYTHON_BINDINGS=ON \
  -D CMAKE_PREFIX_PATH=<path-to-installed-dependencies-if-not-in-default-paths>
```

This plugin is known to build without errors with (as of September 2022)

- [GCC](https://gcc.gnu.org) versions 8 and 9
- [CUDA](https://developer.nvidia.com/cuda-zone) v11.5.1
- [libtorch](https://pytorch.org/) v1.10.2 for CUDA version 10.2 and cxx-11-abi ([download](https://download.pytorch.org/libtorch/cu102/libtorch-cxx11-abi-shared-with-deps-1.10.2%2Bcu102.zip))

*For the ONNX backend:*
- [ONNX](https://github.com/microsoft/onnxruntime) v1.10.0 with CUDA support enabled

*For the Torch backend*
- [TorchScatter](https://github.com/rusty1s/pytorch_scatter) v2.0.9

There were experienced problems with recent GCC 11 versions and CUDA 11.6. A docker image with all dependencies can be found [here](https://github.com/acts-project/machines).

## Running

The Examples of this plugin provide a python-script using the python-bindings to demonstrate the use of the track finding module. The script can be found here:

```bash
./Examples/Scripts/Python/Gnn.py          # tries to run torch
./Examples/Scripts/Python/Gnn.py onnx     # tries to run onnx
./Examples/Scripts/Python/Gnn.py torch    # tries to run torch
```

In order that python can find the `acts.examples` module, set up your `PYTHONPATH` with `source <build>/python/setup.sh`.

## Required files

The track finding module requires some ONNX-files or TorchScript-files that describe the used neural networks. These files are currently not provided within the ACTS repository.

// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsPlugins/Gnn/BoostTrackBuilding.hpp"
#include "ActsPlugins/Gnn/CudaTrackBuilding.hpp"
#include "ActsPlugins/Gnn/GnnPipeline.hpp"
#include "ActsPlugins/Gnn/ModuleMapCuda.hpp"
#include "ActsPlugins/Gnn/OnnxEdgeClassifier.hpp"
#include "ActsPlugins/Gnn/TensorRTEdgeClassifier.hpp"
#include "ActsPlugins/Gnn/TorchEdgeClassifier.hpp"
#include "ActsPlugins/Gnn/TorchMetricLearning.hpp"
#include "ActsPlugins/Gnn/TruthGraphMetricsHook.hpp"
#include "ActsPython/Utilities/Macros.hpp"

#include <boost/preprocessor/if.hpp>
#include <boost/vmd/tuple/size.hpp>
#include <pybind11/functional.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#define ACTS_PYTHON_DECLARE_GNN_STAGE(algorithm, base, mod, ...)            \
  do {                                                                      \
    using namespace Acts;                                                   \
                                                                            \
    using Alg = algorithm;                                                  \
    using Config = Alg::Config;                                             \
    auto alg = py::class_<Alg, base, std::shared_ptr<Alg>>(mod, #algorithm) \
                   .def(py::init([](const Config &c, Logging::Level lvl) {  \
                          return std::make_shared<Alg>(                     \
                              c, getDefaultLogger(#algorithm, lvl));        \
                        }),                                                 \
                        py::arg("config"), py::arg("level"))                \
                   .def_property_readonly("config", &Alg::config);          \
                                                                            \
    auto c = py::class_<Config>(alg, "Config").def(py::init<>());           \
    BOOST_PP_IF(BOOST_VMD_IS_EMPTY(__VA_ARGS__), BOOST_PP_EMPTY(),          \
                ACTS_PYTHON_STRUCT(c, __VA_ARGS__));                        \
  } while (0)

namespace py = pybind11;

using namespace Acts;
using namespace ActsPlugins;
using namespace ActsPython;
using namespace py::literals;

PYBIND11_MODULE(ActsPluginsPythonBindingsGnn, gnn) {
  {
    using C = GraphConstructionBase;
    auto c = py::class_<C, std::shared_ptr<C>>(gnn, "GraphConstructionBase");
  }
  {
    using C = EdgeClassificationBase;
    auto c = py::class_<C, std::shared_ptr<C>>(gnn, "EdgeClassificationBase");
  }
  {
    using C = TrackBuildingBase;
    auto c = py::class_<C, std::shared_ptr<C>>(gnn, "TrackBuildingBase");
  }

  ACTS_PYTHON_DECLARE_GNN_STAGE(BoostTrackBuilding, TrackBuildingBase, gnn);

#ifdef ACTS_GNN_TORCH_BACKEND
  ACTS_PYTHON_DECLARE_GNN_STAGE(TorchMetricLearning, GraphConstructionBase, gnn,
                                modelPath, selectedFeatures, embeddingDim, rVal,
                                knnVal, deviceID);

  ACTS_PYTHON_DECLARE_GNN_STAGE(TorchEdgeClassifier, EdgeClassificationBase,
                                gnn, modelPath, selectedFeatures, cut, nChunks,
                                undirected, deviceID, useEdgeFeatures);
#endif

#ifdef ACTS_GNN_WITH_TENSORRT
  ACTS_PYTHON_DECLARE_GNN_STAGE(TensorRTEdgeClassifier, EdgeClassificationBase,
                                gnn, modelPath, selectedFeatures, cut,
                                numExecutionContexts);
#endif

#ifdef ACTS_GNN_WITH_CUDA
  ACTS_PYTHON_DECLARE_GNN_STAGE(CudaTrackBuilding, TrackBuildingBase, gnn,
                                useOneBlockImplementation, doJunctionRemoval,
                                minCandidateSize);
#endif

#ifdef ACTS_GNN_ONNX_BACKEND
  ACTS_PYTHON_DECLARE_GNN_STAGE(OnnxEdgeClassifier, EdgeClassificationBase, gnn,
                                modelPath, cut);
#endif

#ifdef ACTS_GNN_WITH_MODULEMAP
  ACTS_PYTHON_DECLARE_GNN_STAGE(
      ModuleMapCuda, GraphConstructionBase, gnn, moduleMapPath, rScale,
      phiScale, zScale, etaScale, moreParallel, gpuDevice, gpuBlocks, epsilon);
#endif

  {
    auto cls = py::class_<GnnHook, std::shared_ptr<GnnHook>>(gnn, "GnnHook");
  }

  {
    using Class = TruthGraphMetricsHook;

    auto cls = py::class_<Class, GnnHook, std::shared_ptr<Class>>(
                   gnn, "TruthGraphMetricsHook")
                   .def(py::init([](const std::vector<std::int64_t> &g,
                                    Logging::Level lvl) {
                     return std::make_shared<Class>(
                         g, getDefaultLogger("TruthGraphHook", lvl));
                   }));
  }

  {
    auto cls = py::class_<Device>(gnn, "Device")
                   .def_static("Cpu", &Device::Cpu)
                   .def_static("Cuda", &Device::Cuda, py::arg("index") = 0);
  }

  {
    using Class = GnnPipeline;

    auto cls =
        py::class_<Class, std::shared_ptr<Class>>(gnn, "GnnPipeline")
            .def(py::init(
                     [](std::shared_ptr<GraphConstructionBase> g,
                        std::vector<std::shared_ptr<EdgeClassificationBase>> e,
                        std::shared_ptr<TrackBuildingBase> t,
                        Logging::Level lvl) {
                       return std::make_shared<Class>(
                           g, e, t, getDefaultLogger("MetricLearning", lvl));
                     }),
                 py::arg("graphConstructor"), py::arg("edgeClassifiers"),
                 py::arg("trackBuilder"), py::arg("level"))
            .def("run", &GnnPipeline::run, py::arg("features"),
                 py::arg("moduleIds"), py::arg("spacepoints"),
                 py::arg("device") = Device::Cuda(0),
                 py::arg("hook") = GnnHook{}, py::arg("timing") = nullptr);
  }
}

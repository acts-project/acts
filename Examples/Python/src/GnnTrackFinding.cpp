// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Plugins/Gnn/BoostTrackBuilding.hpp"
#include "Acts/Plugins/Gnn/CudaTrackBuilding.hpp"
#include "Acts/Plugins/Gnn/GnnPipeline.hpp"
#include "Acts/Plugins/Gnn/ModuleMapCuda.hpp"
#include "Acts/Plugins/Gnn/OnnxEdgeClassifier.hpp"
#include "Acts/Plugins/Gnn/TensorRTEdgeClassifier.hpp"
#include "Acts/Plugins/Gnn/TorchEdgeClassifier.hpp"
#include "Acts/Plugins/Gnn/TorchMetricLearning.hpp"
#include "Acts/Plugins/Gnn/TruthGraphMetricsHook.hpp"
#include "ActsExamples/TrackFindingGnn/PrototracksToParameters.hpp"
#include "ActsExamples/TrackFindingGnn/TrackFindingAlgorithmGnn.hpp"
#include "ActsExamples/TrackFindingGnn/TrackFindingFromPrototrackAlgorithm.hpp"
#include "ActsExamples/TrackFindingGnn/TruthGraphBuilder.hpp"
#include "ActsPython/Utilities/Helpers.hpp"
#include "ActsPython/Utilities/Macros.hpp"

#include <memory>

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

using namespace ActsExamples;
using namespace Acts;
using namespace py::literals;

namespace ActsPython {

void addGnnTrackFinding(Context &ctx) {
  auto [m, mex] = ctx.get("main", "examples");

  {
    using C = GraphConstructionBase;
    auto c = py::class_<C, std::shared_ptr<C>>(mex, "GraphConstructionBase");
  }
  {
    using C = EdgeClassificationBase;
    auto c = py::class_<C, std::shared_ptr<C>>(mex, "EdgeClassificationBase");
  }
  {
    using C = TrackBuildingBase;
    auto c = py::class_<C, std::shared_ptr<C>>(mex, "TrackBuildingBase");
  }

  ACTS_PYTHON_DECLARE_GNN_STAGE(BoostTrackBuilding, TrackBuildingBase, mex);

#ifdef ACTS_GNN_TORCH_BACKEND
  ACTS_PYTHON_DECLARE_GNN_STAGE(TorchMetricLearning, GraphConstructionBase, mex,
                                modelPath, selectedFeatures, embeddingDim, rVal,
                                knnVal, deviceID);

  ACTS_PYTHON_DECLARE_GNN_STAGE(TorchEdgeClassifier, EdgeClassificationBase,
                                mex, modelPath, selectedFeatures, cut, nChunks,
                                undirected, deviceID, useEdgeFeatures);
#endif

#ifdef ACTS_GNN_WITH_TENSORRT
  ACTS_PYTHON_DECLARE_GNN_STAGE(TensorRTEdgeClassifier, EdgeClassificationBase,
                                mex, modelPath, selectedFeatures, cut,
                                numExecutionContexts);
#endif

#ifdef ACTS_GNN_WITH_CUDA
  ACTS_PYTHON_DECLARE_GNN_STAGE(CudaTrackBuilding, TrackBuildingBase, mex,
                                useOneBlockImplementation, doJunctionRemoval);
#endif

#ifdef ACTS_GNN_ONNX_BACKEND
  ACTS_PYTHON_DECLARE_GNN_STAGE(OnnxEdgeClassifier, EdgeClassificationBase, mex,
                                modelPath, cut);
#endif

#ifdef ACTS_GNN_WITH_MODULEMAP
  ACTS_PYTHON_DECLARE_GNN_STAGE(
      ModuleMapCuda, GraphConstructionBase, mex, moduleMapPath, rScale,
      phiScale, zScale, etaScale, moreParallel, gpuDevice, gpuBlocks, epsilon);
#endif

  ACTS_PYTHON_DECLARE_ALGORITHM(TruthGraphBuilder, mex, "TruthGraphBuilder",
                                inputSpacePoints, inputSimHits, inputParticles,
                                inputMeasurementSimHitsMap,
                                inputMeasurementParticlesMap, outputGraph,
                                targetMinPT, targetMinSize, uniqueModules);

  {
    auto nodeFeatureEnum =
        py::enum_<TrackFindingAlgorithmGnn::NodeFeature>(mex, "NodeFeature")
            .value("R", TrackFindingAlgorithmGnn::NodeFeature::eR)
            .value("Phi", TrackFindingAlgorithmGnn::NodeFeature::ePhi)
            .value("Z", TrackFindingAlgorithmGnn::NodeFeature::eZ)
            .value("X", TrackFindingAlgorithmGnn::NodeFeature::eX)
            .value("Y", TrackFindingAlgorithmGnn::NodeFeature::eY)
            .value("Eta", TrackFindingAlgorithmGnn::NodeFeature::eEta)
            .value("ClusterX",
                   TrackFindingAlgorithmGnn::NodeFeature::eClusterLoc0)
            .value("ClusterY",
                   TrackFindingAlgorithmGnn::NodeFeature::eClusterLoc1)
            .value("CellCount",
                   TrackFindingAlgorithmGnn::NodeFeature::eCellCount)
            .value("ChargeSum",
                   TrackFindingAlgorithmGnn::NodeFeature::eChargeSum);

    // clang-format off
#define ADD_FEATURE_ENUMS(n) \
  nodeFeatureEnum \
    .value("Cluster" #n "X", TrackFindingAlgorithmGnn::NodeFeature::eCluster##n##X) \
    .value("Cluster" #n "Y", TrackFindingAlgorithmGnn::NodeFeature::eCluster##n##Y) \
    .value("Cluster" #n "Z", TrackFindingAlgorithmGnn::NodeFeature::eCluster##n##Z) \
    .value("Cluster" #n "R", TrackFindingAlgorithmGnn::NodeFeature::eCluster##n##R) \
    .value("Cluster" #n "Phi", TrackFindingAlgorithmGnn::NodeFeature::eCluster##n##Phi) \
    .value("Cluster" #n "Eta", TrackFindingAlgorithmGnn::NodeFeature::eCluster##n##Eta) \
    .value("CellCount" #n, TrackFindingAlgorithmGnn::NodeFeature::eCellCount##n) \
    .value("ChargeSum" #n, TrackFindingAlgorithmGnn::NodeFeature::eChargeSum##n) \
    .value("LocEta" #n, TrackFindingAlgorithmGnn::NodeFeature::eLocEta##n) \
    .value("LocPhi" #n, TrackFindingAlgorithmGnn::NodeFeature::eLocPhi##n) \
    .value("LocDir0" #n, TrackFindingAlgorithmGnn::NodeFeature::eLocDir0##n) \
    .value("LocDir1" #n, TrackFindingAlgorithmGnn::NodeFeature::eLocDir1##n) \
    .value("LocDir2" #n, TrackFindingAlgorithmGnn::NodeFeature::eLocDir2##n) \
    .value("LengthDir0" #n, TrackFindingAlgorithmGnn::NodeFeature::eLengthDir0##n) \
    .value("LengthDir1" #n, TrackFindingAlgorithmGnn::NodeFeature::eLengthDir1##n) \
    .value("LengthDir2" #n, TrackFindingAlgorithmGnn::NodeFeature::eLengthDir2##n) \
    .value("GlobEta" #n, TrackFindingAlgorithmGnn::NodeFeature::eGlobEta##n) \
    .value("GlobPhi" #n, TrackFindingAlgorithmGnn::NodeFeature::eGlobPhi##n) \
    .value("EtaAngle" #n, TrackFindingAlgorithmGnn::NodeFeature::eEtaAngle##n) \
    .value("PhiAngle" #n, TrackFindingAlgorithmGnn::NodeFeature::ePhiAngle##n)
    // clang-format on

    ADD_FEATURE_ENUMS(1);
    ADD_FEATURE_ENUMS(2);

#undef ADD_FEATURE_ENUMS
  }

  ACTS_PYTHON_DECLARE_ALGORITHM(
      TrackFindingAlgorithmGnn, mex, "TrackFindingAlgorithmGnn",
      inputSpacePoints, inputClusters, inputTruthGraph, outputProtoTracks,
      outputGraph, graphConstructor, edgeClassifiers, trackBuilder,
      nodeFeatures, featureScales, minMeasurementsPerTrack, geometryIdMap);

  { auto cls = py::class_<GnnHook, std::shared_ptr<GnnHook>>(mex, "GnnHook"); }

  {
    using Class = TruthGraphMetricsHook;

    auto cls = py::class_<Class, GnnHook, std::shared_ptr<Class>>(
                   mex, "TruthGraphMetricsHook")
                   .def(py::init([](const std::vector<std::int64_t> &g,
                                    Logging::Level lvl) {
                     return std::make_shared<Class>(
                         g, getDefaultLogger("TruthGraphHook", lvl));
                   }));
  }

  {
    auto cls = py::class_<Device>(mex, "Device")
                   .def_static("Cpu", &Device::Cpu)
                   .def_static("Cuda", &Device::Cuda, py::arg("index") = 0);
  }

  {
    using Class = GnnPipeline;

    auto cls =
        py::class_<Class, std::shared_ptr<Class>>(mex, "GnnPipeline")
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

  ACTS_PYTHON_DECLARE_ALGORITHM(
      PrototracksToParameters, mex, "PrototracksToParameters", inputProtoTracks,
      inputSpacePoints, outputSeeds, outputParameters, outputProtoTracks,
      geometry, magneticField, buildTightSeeds);

  ACTS_PYTHON_DECLARE_ALGORITHM(
      TrackFindingFromPrototrackAlgorithm, mex,
      "TrackFindingFromPrototrackAlgorithm", inputProtoTracks,
      inputMeasurements, inputInitialTrackParameters, outputTracks,
      measurementSelectorCfg, trackingGeometry, magneticField, findTracks, tag);
}

}  // namespace ActsPython

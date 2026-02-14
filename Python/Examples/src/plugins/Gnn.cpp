// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/TrackFindingGnn/ProtoTracksToParameters.hpp"
#include "ActsExamples/TrackFindingGnn/TrackFindingAlgorithmGnn.hpp"
#include "ActsExamples/TrackFindingGnn/TrackFindingFromProtoTracksAlgorithm.hpp"
#include "ActsExamples/TrackFindingGnn/TruthGraphBuilder.hpp"
#include "ActsPlugins/Gnn/BoostTrackBuilding.hpp"
#include "ActsPlugins/Gnn/CudaTrackBuilding.hpp"
#include "ActsPlugins/Gnn/GnnPipeline.hpp"
#include "ActsPlugins/Gnn/ModuleMapCuda.hpp"
#include "ActsPlugins/Gnn/OnnxEdgeClassifier.hpp"
#include "ActsPlugins/Gnn/TensorRTEdgeClassifier.hpp"
#include "ActsPlugins/Gnn/TorchEdgeClassifier.hpp"
#include "ActsPlugins/Gnn/TorchMetricLearning.hpp"
#include "ActsPlugins/Gnn/TruthGraphMetricsHook.hpp"
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

using namespace Acts;
using namespace ActsPlugins;
using namespace ActsExamples;
using namespace ActsPython;
using namespace py::literals;

PYBIND11_MODULE(ActsExamplesPythonBindingsGnn, gnn) {
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
                                nugnnecutionContexts);
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

  ACTS_PYTHON_DECLARE_ALGORITHM(TruthGraphBuilder, gnn, "TruthGraphBuilder",
                                inputSpacePoints, inputSimHits, inputParticles,
                                inputMeasurementSimHitsMap,
                                inputMeasurementParticlesMap, outputGraph,
                                targetMinPT, targetMinSize, uniqueModules);

  {
    auto nodeFeatureEnum =
        py::enum_<TrackFindingAlgorithmGnn::NodeFeature>(gnn, "NodeFeature")
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
      TrackFindingAlgorithmGnn, gnn, "TrackFindingAlgorithmGnn",
      inputSpacePoints, inputClusters, inputTruthGraph, outputProtoTracks,
      outputGraph, graphConstructor, edgeClassifiers, trackBuilder,
      nodeFeatures, featureScales, minMeasurementsPerTrack, geometryIdMap);

  { auto cls = py::class_<GnnHook, std::shared_ptr<GnnHook>>(gnn, "GnnHook"); }

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

  ACTS_PYTHON_DECLARE_ALGORITHM(
      ProtoTracksToParameters, gnn, "ProtoTracksToParameters", inputProtoTracks,
      inputSpacePoints, outputSeeds, outputParameters, outputProtoTracks,
      geometry, magneticField, buildTightSeeds);

  ACTS_PYTHON_DECLARE_ALGORITHM(
      TrackFindingFromProtoTracksAlgorithm, gnn,
      "TrackFindingFromProtoTracksAlgorithm", inputProtoTracks,
      inputMeasurements, inputInitialTrackParameters, outputTracks,
      measurementSelectorCfg, trackingGeometry, magneticField, findTracks, tag);
}

// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Plugins/ExaTrkX/BoostTrackBuilding.hpp"
#include "Acts/Plugins/ExaTrkX/CudaTrackBuilding.hpp"
#include "Acts/Plugins/ExaTrkX/DummyEdgeClassifier.hpp"
#include "Acts/Plugins/ExaTrkX/ExaTrkXPipeline.hpp"
#include "Acts/Plugins/ExaTrkX/FullyConnectedGraphConstructor.hpp"
#include "Acts/Plugins/ExaTrkX/ModuleMapCpp.hpp"
#include "Acts/Plugins/ExaTrkX/ModuleMapCuda.hpp"
#include "Acts/Plugins/ExaTrkX/OnnxEdgeClassifier.hpp"
#include "Acts/Plugins/ExaTrkX/TensorRTEdgeClassifier.hpp"
#include "Acts/Plugins/ExaTrkX/TorchEdgeClassifier.hpp"
#include "Acts/Plugins/ExaTrkX/TorchEdgeClassifierAOT.hpp"
#include "Acts/Plugins/ExaTrkX/TorchMetricLearning.hpp"
#include "Acts/Plugins/ExaTrkX/TorchTruthGraphMetricsHook.hpp"
#include "Acts/Plugins/Python/Utilities.hpp"
#include "ActsExamples/EventData/GeometryContainers.hpp"
#include "ActsExamples/TrackFindingExaTrkX/PrototracksToParameters.hpp"
#include "ActsExamples/TrackFindingExaTrkX/TrackFindingAlgorithmExaTrkX.hpp"
#include "ActsExamples/TrackFindingExaTrkX/TrackFindingFromPrototrackAlgorithm.hpp"
#include "ActsExamples/TrackFindingExaTrkX/TruthGraphBuilder.hpp"

#include <memory>

#include <boost/algorithm/string.hpp>
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
    ACTS_PYTHON_STRUCT(c, __VA_ARGS__);                                     \
  } while (0)

namespace py = pybind11;

using namespace ActsExamples;
using namespace Acts;
using namespace py::literals;

namespace Acts::Python {

void addExaTrkXTrackFinding(Context &ctx) {
  auto [m, mex] = ctx.get("main", "examples");

  {
    using C = Acts::GraphConstructionBase;
    auto c = py::class_<C, std::shared_ptr<C>>(mex, "GraphConstructionBase");
  }
  {
    using C = Acts::EdgeClassificationBase;
    auto c = py::class_<C, std::shared_ptr<C>>(mex, "EdgeClassificationBase");
  }
  {
    using C = Acts::TrackBuildingBase;
    auto c = py::class_<C, std::shared_ptr<C>>(mex, "TrackBuildingBase");
  }

  ACTS_PYTHON_DECLARE_GNN_STAGE(DummyEdgeClassifier, EdgeClassificationBase,
                                mex, keepAll);

  ACTS_PYTHON_DECLARE_GNN_STAGE(BoostTrackBuilding, TrackBuildingBase, mex,
                                doJunctionRemoval, doWalkthrough,
                                walkthroughLowCut, walkthroughHighCut);

#ifdef ACTS_EXATRKX_TORCH_BACKEND
  ACTS_PYTHON_DECLARE_GNN_STAGE(
      FullyConnectedGraphConstructor, GraphConstructionBase, mex, maxGraphSize,
      maxOutEdges, maxDeltaR, maxDeltaZ, rScale, phiScale, zScale, etaScale,
      rOffset, phiOffset, zOffset, etaOffset);

  ACTS_PYTHON_DECLARE_GNN_STAGE(TorchMetricLearning, GraphConstructionBase, mex,
                                modelPath, selectedFeatures, embeddingDim, rVal,
                                knnVal, deviceID);

  ACTS_PYTHON_DECLARE_GNN_STAGE(TorchEdgeClassifier, EdgeClassificationBase,
                                mex, modelPath, selectedFeatures, cut, nChunks,
                                undirected, deviceID, useEdgeFeatures);
#endif

#ifdef ACTS_EXATRKX_WITH_TENSORRT
  ACTS_PYTHON_DECLARE_GNN_STAGE(TensorRTEdgeClassifier, EdgeClassificationBase,
                                mex, modelPath, selectedFeatures, cut,
                                numExecutionContexts);
#endif

#ifdef ACTS_EXATRKX_WITH_CUDA
  ACTS_PYTHON_DECLARE_GNN_STAGE(CudaTrackBuilding, TrackBuildingBase, mex,
                                useOneBlockImplementation, doJunctionRemoval);
#endif

#ifdef ACTS_EXATRKX_ONNX_BACKEND
  ACTS_PYTHON_DECLARE_GNN_STAGE(OnnxEdgeClassifier, EdgeClassificationBase, mex,
                                modelPath, cut);
#endif

#ifdef ACTS_EXATRKX_WITH_MODULEMAP
  ACTS_PYTHON_DECLARE_GNN_STAGE(ModuleMapCuda, GraphConstructionBase, mex,
                                moduleMapPath, rScale, phiScale, zScale,
                                gpuDevice, gpuBlocks, moreParallel);
#endif

  ACTS_PYTHON_DECLARE_ALGORITHM(
      ActsExamples::TruthGraphBuilder, mex, "TruthGraphBuilder",
      inputSpacePoints, inputSimHits, inputParticles,
      inputMeasurementSimHitsMap, inputMeasurementParticlesMap, outputGraph,
      targetMinPT, targetMinSize, uniqueModules);

  {
    auto nodeFeatureEnum =
        py::enum_<TrackFindingAlgorithmExaTrkX::NodeFeature>(mex, "NodeFeature")
            .value("R", TrackFindingAlgorithmExaTrkX::NodeFeature::eR)
            .value("Phi", TrackFindingAlgorithmExaTrkX::NodeFeature::ePhi)
            .value("Z", TrackFindingAlgorithmExaTrkX::NodeFeature::eZ)
            .value("X", TrackFindingAlgorithmExaTrkX::NodeFeature::eX)
            .value("Y", TrackFindingAlgorithmExaTrkX::NodeFeature::eY)
            .value("Eta", TrackFindingAlgorithmExaTrkX::NodeFeature::eEta)
            .value("ClusterX",
                   TrackFindingAlgorithmExaTrkX::NodeFeature::eClusterLoc0)
            .value("ClusterY",
                   TrackFindingAlgorithmExaTrkX::NodeFeature::eClusterLoc1)
            .value("CellCount",
                   TrackFindingAlgorithmExaTrkX::NodeFeature::eCellCount)
            .value("ChargeSum",
                   TrackFindingAlgorithmExaTrkX::NodeFeature::eChargeSum);

    // clang-format off
#define ADD_FEATURE_ENUMS(n) \
  nodeFeatureEnum \
    .value("Cluster" #n "X", TrackFindingAlgorithmExaTrkX::NodeFeature::eCluster##n##X) \
    .value("Cluster" #n "Y", TrackFindingAlgorithmExaTrkX::NodeFeature::eCluster##n##Y) \
    .value("Cluster" #n "Z", TrackFindingAlgorithmExaTrkX::NodeFeature::eCluster##n##Z) \
    .value("Cluster" #n "R", TrackFindingAlgorithmExaTrkX::NodeFeature::eCluster##n##R) \
    .value("Cluster" #n "Phi", TrackFindingAlgorithmExaTrkX::NodeFeature::eCluster##n##Phi) \
    .value("Cluster" #n "Eta", TrackFindingAlgorithmExaTrkX::NodeFeature::eCluster##n##Eta) \
    .value("CellCount" #n, TrackFindingAlgorithmExaTrkX::NodeFeature::eCellCount##n) \
    .value("ChargeSum" #n, TrackFindingAlgorithmExaTrkX::NodeFeature::eChargeSum##n) \
    .value("LocEta" #n, TrackFindingAlgorithmExaTrkX::NodeFeature::eLocEta##n) \
    .value("LocPhi" #n, TrackFindingAlgorithmExaTrkX::NodeFeature::eLocPhi##n) \
    .value("LocDir0" #n, TrackFindingAlgorithmExaTrkX::NodeFeature::eLocDir0##n) \
    .value("LocDir1" #n, TrackFindingAlgorithmExaTrkX::NodeFeature::eLocDir1##n) \
    .value("LocDir2" #n, TrackFindingAlgorithmExaTrkX::NodeFeature::eLocDir2##n) \
    .value("LengthDir0" #n, TrackFindingAlgorithmExaTrkX::NodeFeature::eLengthDir0##n) \
    .value("LengthDir1" #n, TrackFindingAlgorithmExaTrkX::NodeFeature::eLengthDir1##n) \
    .value("LengthDir2" #n, TrackFindingAlgorithmExaTrkX::NodeFeature::eLengthDir2##n) \
    .value("GlobEta" #n, TrackFindingAlgorithmExaTrkX::NodeFeature::eGlobEta##n) \
    .value("GlobPhi" #n, TrackFindingAlgorithmExaTrkX::NodeFeature::eGlobPhi##n) \
    .value("EtaAngle" #n, TrackFindingAlgorithmExaTrkX::NodeFeature::eEtaAngle##n) \
    .value("PhiAngle" #n, TrackFindingAlgorithmExaTrkX::NodeFeature::ePhiAngle##n)
    // clang-format on

    ADD_FEATURE_ENUMS(1);
    ADD_FEATURE_ENUMS(2);

#undef ADD_FEATURE_ENUMS
  }

  ACTS_PYTHON_DECLARE_ALGORITHM(
      ActsExamples::TrackFindingAlgorithmExaTrkX, mex,
      "TrackFindingAlgorithmExaTrkX", inputSpacePoints, inputClusters,
      inputTruthGraph, outputProtoTracks, outputGraph, graphConstructor,
      edgeClassifiers, trackBuilder, nodeFeatures, featureScales,
      minMeasurementsPerTrack, geometryIdMap);

  {
    auto cls =
        py::class_<Acts::ExaTrkXHook, std::shared_ptr<Acts::ExaTrkXHook>>(
            mex, "ExaTrkXHook");
  }

  {
    using Class = Acts::TorchTruthGraphMetricsHook;

    auto cls = py::class_<Class, Acts::ExaTrkXHook, std::shared_ptr<Class>>(
                   mex, "TorchTruthGraphMetricsHook")
                   .def(py::init([](const std::vector<std::int64_t> &g,
                                    Logging::Level lvl) {
                     return std::make_shared<Class>(
                         g, getDefaultLogger("PipelineHook", lvl));
                   }));
  }

  {
    using Class = Acts::ExaTrkXPipeline;

    auto cls =
        py::class_<Class, std::shared_ptr<Class>>(mex, "ExaTrkXPipeline")
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
            .def("run", &ExaTrkXPipeline::run, py::arg("features"),
                 py::arg("moduleIds"), py::arg("spacepoints"),
                 py::arg("hook") = Acts::ExaTrkXHook{},
                 py::arg("timing") = nullptr);
  }

  ACTS_PYTHON_DECLARE_ALGORITHM(
      ActsExamples::PrototracksToParameters, mex, "PrototracksToParameters",
      inputProtoTracks, inputSpacePoints, outputSeeds, outputParameters,
      outputProtoTracks, geometry, magneticField, buildTightSeeds,
      minSpacepointDist, stripVolumes, initialSigmas);

  ACTS_PYTHON_DECLARE_ALGORITHM(
      ActsExamples::TrackFindingFromPrototrackAlgorithm, mex,
      "TrackFindingFromPrototrackAlgorithm", inputProtoTracks,
      inputMeasurements, inputInitialTrackParameters, outputTracks,
      measurementSelectorCfg, trackingGeometry, magneticField, findTracks, tag,
      onlyPrototrackMeasurements, pickSeed);

  py::class_<GeometryIdMapActsAthena, std::shared_ptr<GeometryIdMapActsAthena>>(
      mex, "GeometryIdMapActsAthena")
      .def(py::init<>())
      .def_static(
          "fromCsv",
          [](const std::string &filename) {
            std::ifstream inputfile(filename);
            if (!inputfile) {
              throw std::invalid_argument("Could not open '" + filename + "'");
            }

            auto map = std::make_shared<GeometryIdMapActsAthena>();

            // Read the header first and don't process it
            std::string line;
            std::getline(inputfile, line);

            // Read the tokens in the vector
            std::vector<std::string> tokens;
            for (; std::getline(inputfile, line);) {
              boost::split(tokens, line, boost::is_any_of(","));
              auto athenaId = std::stoul(tokens.at(0));
              auto actsId = std::stoul(tokens.at(1));
              map->left.insert({athenaId, Acts::GeometryIdentifier{actsId}});
              tokens.clear();
            }

            return map;
          },
          "filename"_a);
}

}  // namespace Acts::Python

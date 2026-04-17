// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/TrackFindingGnn/TrackFindingAlgorithmGnn.hpp"
#include "ActsExamples/TrackFindingGnn/TrackFindingFromProtoTracksAlgorithm.hpp"
#include "ActsExamples/TrackFindingGnn/TruthGraphBuilder.hpp"
#include "ActsPython/Utilities/Macros.hpp"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;

using namespace Acts;
using namespace ActsExamples;
using namespace ActsPython;
using namespace py::literals;

PYBIND11_MODULE(ActsExamplesPythonBindingsGnn, gnn) {
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
      nodeFeatures, featureScales, minMeasurementsPerTrack, geometryIdMap,
      device);

  ACTS_PYTHON_DECLARE_ALGORITHM(
      TrackFindingFromProtoTracksAlgorithm, gnn,
      "TrackFindingFromProtoTracksAlgorithm", inputProtoTracks,
      inputMeasurements, inputInitialTrackParameters, outputTracks,
      measurementSelectorCfg, trackingGeometry, magneticField, findTracks, tag);
}

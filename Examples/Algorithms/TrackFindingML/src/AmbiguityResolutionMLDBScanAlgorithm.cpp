// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/TrackFindingML/AmbiguityResolutionMLDBScanAlgorithm.hpp"

#include "Acts/Plugins/Mlpack/AmbiguityDBScanClustering.hpp"
#include "ActsExamples/Framework/ProcessCode.hpp"
#include "ActsExamples/Framework/WhiteBoard.hpp"

#include <iterator>
#include <map>

ActsExamples::AmbiguityResolutionMLDBScanAlgorithm::
    AmbiguityResolutionMLDBScanAlgorithm(
        ActsExamples::AmbiguityResolutionMLDBScanAlgorithm::Config cfg,
        Acts::Logging::Level lvl)
    : ActsExamples::AmbiguityResolutionML(
          "AmbiguityResolutionMLDBScanAlgorithm", lvl),
      m_cfg(std::move(cfg)),
      m_duplicateClassifier(m_cfg.inputDuplicateNN.c_str()) {
  if (m_cfg.inputTracks.empty()) {
    throw std::invalid_argument("Missing trajectories input collection");
  }
  if (m_cfg.outputTracks.empty()) {
    throw std::invalid_argument("Missing trajectories output collection");
  }
  m_inputTracks.initialize(m_cfg.inputTracks);
  m_outputTracks.initialize(m_cfg.outputTracks);
}

ActsExamples::ProcessCode
ActsExamples::AmbiguityResolutionMLDBScanAlgorithm::execute(
    const AlgorithmContext& ctx) const {
  // Read input data
  const auto& tracks = m_inputTracks(ctx);
  // Associate measurement to their respective tracks
  std::multimap<int, std::pair<int, std::vector<int>>> trackMap =
      mapTrackHits(tracks, m_cfg.nMeasurementsMin);
  // Cluster the tracks using DBscan
  auto cluster = Acts::dbscanTrackClustering(
      trackMap, tracks, m_cfg.epsilonDBScan, m_cfg.minPointsDBScan);
  // Select the ID of the track we want to keep
  std::vector<int> goodTracks =
      m_duplicateClassifier.solveAmbuguity(cluster, tracks);
  // Prepare the output track collection from the IDs
  auto outputTracks = prepareOutputTrack(tracks, goodTracks);
  m_outputTracks(ctx, std::move(outputTracks));

  return ActsExamples::ProcessCode::SUCCESS;
}

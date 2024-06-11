// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Units.hpp"
#include "Acts/EventData/TrackContainer.hpp"
#include "Acts/Utilities/Delegate.hpp"
#include "Acts/Utilities/Logger.hpp"

#include <cstddef>
#include <map>
#include <memory>
#include <string>
#include <tuple>
#include <vector>


namespace Acts {

/// Generic implementation of the machine learning ambiguity resolution
/// Contains method for data preparations
template <typename AmbiguityNetwork>
class AmbiguityResolutionML {
 public:
  struct Config {
    /// Path to the ONNX model for the duplicate neural network
    std::string inputDuplicateNN;
    /// Minimum number of measurement to form a track.
    std::size_t nMeasurementsMin = 7;
  };
  /// Construct the ambiguity resolution algorithm.
  ///
  /// @param name name of the algorithm
  /// @param lvl is the logging level
  AmbiguityResolutionML(const Config& cfg,
                        std::unique_ptr<const Logger> logger = getDefaultLogger(
                            "AmbiguityResolutionML", Logging::INFO))
      : m_cfg{cfg},
        m_duplicateClassifier(m_cfg.inputDuplicateNN.c_str()),
        m_logger{std::move(logger)} {}

  /// Associated measurements ID to Tracks ID
  ///
  /// @param tracks is the input track container
  /// @param nMeasurementsMin minimum number of measurement per track
  /// @return an ordered list containing pairs of track ID and associated measurement ID
  template <typename track_container_t, typename traj_t,
            template <typename> class holder_t, typename source_link_hash_t,
            typename source_link_equality_t>
  std::multimap<int, std::pair<std::size_t, std::vector<std::size_t>>>
  mapTrackHits(
      const TrackContainer<track_container_t, traj_t, holder_t>& tracks,
      source_link_hash_t&& sourceLinkHash,
      source_link_equality_t&& sourceLinkEquality) const {
    auto measurementIndexMap =
        std::unordered_map<SourceLink, std::size_t, source_link_hash_t,
                           source_link_equality_t>(0, sourceLinkHash,
                                                   sourceLinkEquality);

    std::multimap<int, std::pair<std::size_t, std::vector<std::size_t>>>
        trackMap;
    std::size_t trackIndex = 0;
    // Loop over all the trajectories in the events
    for (const auto& track : tracks) {
      // Kick out tracks that do not fulfill our initial requirements
      if (track.nMeasurements() < m_cfg.nMeasurementsMin) {
        continue;
      }
      std::vector<std::size_t> measurements;
      for (auto ts : track.trackStatesReversed()) {
        if (ts.typeFlags().test(Acts::TrackStateFlag::MeasurementFlag)) {
          SourceLink sourceLink = ts.getUncalibratedSourceLink();
          // assign a new measurement index if the source link was not seen yet
          auto emplace = measurementIndexMap.try_emplace(
              sourceLink, measurementIndexMap.size());
          measurements.push_back(emplace.first->second);
        }
      }
      trackMap.emplace(track.nMeasurements(),
                       std::make_pair(trackIndex, measurements));
      ++trackIndex;
    }
    return trackMap;
  }

  /// Select the track associated with each cluster
  ///
  /// @param clusters is a map of clusters, each cluster correspond to a vector of track ID
  /// @param tracks is the input track container
  /// @return a vector of trackID corresponding tho the good tracks
  template <typename track_container_t, typename traj_t,
            template <typename> class holder_t>
  std::vector<std::size_t> solveAmbiguity(
      std::unordered_map<std::size_t, std::vector<std::size_t>>& clusters,
      const Acts::TrackContainer<track_container_t, traj_t, holder_t>& tracks)
      const {
    std::vector<std::vector<float>> outputTensor =
        m_duplicateClassifier.inferScores(clusters, tracks);
    std::vector<std::size_t> goodTracks =
        m_duplicateClassifier.trackSelection(clusters, outputTensor);

    return goodTracks;
  }

 private:
  // Configuration
  Config m_cfg;

  // Onnx model for track scoring
  AmbiguityNetwork m_duplicateClassifier;

  /// Logging instance
  std::unique_ptr<const Logger> m_logger = nullptr;

  /// Private access to logging instance
  const Logger& logger() const;
};

}  // namespace Acts

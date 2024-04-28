// This file is part of the Acts project.
//
// Copyright (C) 2024 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Units.hpp"
#include "Acts/EventData/TrackContainer.hpp"
#include "Acts/Utilities/Delegate.hpp"
#include "Acts/Utilities/Logger.hpp"

#include <map>
#include <string>
#include <vector>

#include <boost/container/flat_map.hpp>
#include <boost/container/flat_set.hpp>
namespace ScoreBasedAmbiguitySolver {
// std::functions defined, to be used in the optional cuts.
template <typename track_container_t, typename trajectory_t,
          template <typename> class holder_t>
using OptionalFilter = std::function<bool(
    const Acts::TrackProxy<track_container_t, trajectory_t, holder_t, true>&)>;
template <typename track_container_t, typename trajectory_t,
          template <typename> class holder_t>
using OptionalScoreModifier = std::function<void(
    const Acts::TrackProxy<track_container_t, trajectory_t, holder_t, true>&,
    double&)>;

}  // namespace ScoreBasedAmbiguitySolver

namespace Acts {

/// Generic implementation of the score based ambiguity resolution.
/// The alhorithm is based on the following steps:
/// 1) Compute the initial state of the tracks
/// 2) Compute the score of each track
/// 3) Removes hits that are not good enough for each track
/// 4) Remove tracks that have a score below a certain threshold or not have
/// enough hits
/// 5) Remove tracks that are not good enough based on cuts Contains method for
/// data preparations
class ScoreBasedAmbiguityResolution {
 public:
  using measurementTuple = std::tuple<std::size_t, std::size_t, bool>;
  /// @brief Detector configuration struct : contains the configuration for each detector
  ///
  /// The configuration can be saved in a json file and loaded from there.
  ///
  struct DetectorConfig {
    int hitsScoreWeight = 0;
    int holesScoreWeight = 0;
    int outliersScoreWeight = 0;
    int otherScoreWeight = 0;

    std::size_t minHits;
    std::size_t maxHits;
    std::size_t maxHoles;
    std::size_t maxOutliers;
    std::size_t maxSharedHits;

    /// if true, the shared hits are considered as bad hits for this detector
    bool sharedHitsFlag;

    std::size_t detectorId;

    /// a list of values from  0 to 1, the higher number of hits, higher value
    /// in the list is multiplied to ambuiguity score applied only if
    /// useAmbiguityFunction is true
    std::vector<double> factorHits;

    /// a list of values from  0 to 1, the higher number of holes, lower value
    /// in the list is multiplied to ambuiguity score applied only if
    /// useAmbiguityFunction is true
    std::vector<double> factorHoles;
  };

  /// @brief  TrackFeatures struct : contains the features that are counted for each track.
  ///
  /// The trackFeatures is used to compute the score of each track
  struct TrackFeatures {
    std::size_t nHits = 0;
    std::size_t nHoles = 0;
    std::size_t nOutliers = 0;
    std::size_t nSharedHits = 0;
  };

  /// @brief Configuration struct : contains the configuration for the ambiguity resolution.
  struct Config {
    std::map<std::size_t, std::size_t> volumeMap = {{0, 0}};
    std::map<std::size_t, DetectorConfig> detectorMap;
    /// minimum score for any track
    double minScore = 0;
    /// minimum score for shared tracks
    double minScoreSharedTracks = 0;
    /// maximum number of shared tracks per measurement
    std::size_t maxSharedTracksPerMeasurement = 10;
    /// maximum number of shared hit per track
    std::size_t maxShared = 5;

    double pTMin = 0 * UnitConstants::GeV;
    double pTMax = 1e5 * UnitConstants::GeV;

    double phiMin = -M_PI * UnitConstants::rad;
    double phiMax = M_PI * UnitConstants::rad;

    double etaMin = -5;
    double etaMax = 5;

    // if true, the ambiguity score is computed based on a different function.
    bool useAmbiguityFunction = false;
  };

  /// @brief Optional_cuts struct : contains the optional cuts to be applied.
  ///
  /// The optional cuts,weights and score are used to remove tracks that are not
  /// good enough, based on some criteria. Users are free to add their own cuts
  /// with the help of this struct.
  template <typename track_container_t, typename traj_t,
            template <typename> class holder_t>
  struct Optional_cuts {
    using OptionalFilter = std::function<bool(
        const Acts::TrackProxy<track_container_t, traj_t, holder_t, true>&)>;

    using OptionalScoreModifier = std::function<void(
        const Acts::TrackProxy<track_container_t, traj_t, holder_t, true>&,
        double&)>;
    std::vector<OptionalFilter> cuts = {};
    std::vector<OptionalScoreModifier> weights = {};
    std::vector<OptionalScoreModifier> ambiscores = {};
    // applied only if useAmbiguityFunction is true
  };
  ScoreBasedAmbiguityResolution(
      const Config& cfg,
      std::unique_ptr<const Logger> logger =
          getDefaultLogger("ScoreBasedAmbiguityResolution", Logging::INFO))
      : m_cfg{cfg}, m_logger{std::move(logger)} {}

  /// Compute the initial state of the tracks.
  ///
  /// @param tracks is the input track container
  /// @param sourceLinkHash is the  source links
  /// @param sourceLinkEquality is the equality function for the source links
  /// @param trackFeaturesMaps is the trackFeatures map from detector ID to trackFeatures
  /// @return a vector of the initial state of the tracks
  template <typename track_container_t, typename traj_t,
            template <typename> class holder_t, typename source_link_hash_t,
            typename source_link_equality_t>
  std::vector<std::vector<measurementTuple>> computeInitialState(
      const TrackContainer<track_container_t, traj_t, holder_t>& tracks,
      source_link_hash_t&& sourceLinkHash,
      source_link_equality_t&& sourceLinkEquality,
      std::vector<std::map<std::size_t, TrackFeatures>>& trackFeaturesMaps)
      const;

  /// Compute the score of each track.
  ///
  /// @param tracks is the input track container
  /// @param trackFeaturesMaps is the trackFeatures map from detector ID to trackFeatures
  /// @param optionalCuts is the user defined optional cuts to be applied.
  /// @return a vector of scores for each track
  template <typename track_container_t, typename traj_t,
            template <typename> class holder_t>
  std::vector<double> simpleScore(
      const TrackContainer<track_container_t, traj_t, holder_t>& tracks,
      const std::vector<std::map<std::size_t, TrackFeatures>>&
          trackFeaturesMaps,
      const Optional_cuts<track_container_t, traj_t, holder_t>& optionalCuts =
          {}) const;

  /// Remove hits that are not good enough for each track and removes tracks
  /// that have a score below a certain threshold or not enough hits.
  ///
  /// @brief Remove tracks that are not good enough based on cuts
  /// @param trackScore is the score of each track
  /// @param trackFeaturesMaps is the trackFeatures map for each track
  /// @param measurementsPerTrack is the list of measurements for each track
  /// @return a vector of IDs of the tracks we want to keep
  std::vector<bool> getCleanedOutTracks(
      const std::vector<double>& trackScore,
      const std::vector<std::map<std::size_t, TrackFeatures>>&
          trackFeaturesMaps,
      const std::vector<std::vector<measurementTuple>>& measurementsPerTrack)
      const;

  /// Remove tracks that are bad based on cuts and weighted scores.
  ///
  /// @brief Remove tracks that are not good enough
  /// @param tracks is the input track container
  /// @param measurementsPerTrack is the list of measurements for each track
  /// @param trackFeaturesMaps is the map of detector id to trackFeatures for each track
  /// @param optionalCuts is the optional cuts to be applied
  /// @return a vector of IDs of the tracks we want to keep
  template <typename track_container_t, typename traj_t,
            template <typename> class holder_t>
  std::vector<int> solveAmbiguity(
      const TrackContainer<track_container_t, traj_t, holder_t>& tracks,
      const std::vector<std::vector<measurementTuple>>& measurementsPerTrack,
      const std::vector<std::map<std::size_t, TrackFeatures>>&
          trackFeaturesMaps,
      const Optional_cuts<track_container_t, traj_t, holder_t>& optionalCuts =
          {}) const;

 private:
  Config m_cfg;

  /// Logging instance
  std::unique_ptr<const Logger> m_logger;

  /// Private access to logging instance
  const Logger& logger() const { return *m_logger; }
};

}  // namespace Acts

#include "Acts/AmbiguityResolution/ScoreBasedAmbiguityResolution.ipp"

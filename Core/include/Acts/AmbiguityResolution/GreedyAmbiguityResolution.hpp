// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/EventData/MultiTrajectoryHelpers.hpp"
#include "Acts/EventData/SourceLink.hpp"
#include "Acts/EventData/TrackContainer.hpp"
#include "Acts/Utilities/Delegate.hpp"
#include "Acts/Utilities/Logger.hpp"

#include <memory>

#include <boost/container/flat_map.hpp>
#include <boost/container/flat_set.hpp>

namespace Acts {

/// Evicts tracks that seem to be duplicated.
///
/// The implementation works as follows:
///  1) Calculate shared hits per track.
///  2) If the maximum shared hits criteria is met, we are done.
///     This is the configurable amount of shared hits we are ok with
///     in our experiment.
///  3) Else, remove the track with the highest relative shared hits (i.e.
///     shared hits / hits).
///  4) Back to square 1.
class GreedyAmbiguityResolution {
 public:
  struct Config {
    /// Maximum amount of shared hits per track.
    std::uint32_t maximumSharedHits = 1;
    /// Maximum number of iterations
    std::uint32_t maximumIterations = 1000;

    /// Minumum number of measurement to form a track.
    size_t nMeasurementsMin = 7;

    /// SourceLink hash.
    Delegate<std::size_t(const SourceLink&)> sourceLinkHash;
    /// SourceLink equality delegate.
    Delegate<bool(const SourceLink&, const SourceLink&)> sourceLinkEquality;
  };

  struct State {
    std::size_t numberOfTracks{};

    std::vector<int> trackTips;
    std::vector<float> trackChi2;
    std::vector<std::vector<std::size_t>> measurementsPerTrack;

    boost::container::flat_map<std::size_t,
                               boost::container::flat_set<std::size_t>>
        tracksPerMeasurement;
    std::vector<std::size_t> sharedMeasurementsPerTrack;

    boost::container::flat_set<std::size_t> selectedTracks;
  };

  GreedyAmbiguityResolution(const Config& cfg,
                            std::unique_ptr<const Logger> logger =
                                getDefaultLogger("GreedyAmbiguityResolution",
                                                 Logging::INFO))
      : m_cfg{cfg}, m_logger{std::move(logger)} {}

  template <typename track_container_t, typename traj_t,
            template <typename> class holder_t>
  void computeInitialState(
      const TrackContainer<track_container_t, traj_t, holder_t>& tracks,
      State& state) const;

  void resolve(State& state) const;

 private:
  Config m_cfg;

  /// Logging instance
  std::unique_ptr<const Logger> m_logger;

  /// Private access to logging instance
  const Logger& logger() const { return *m_logger; }

  struct MapableSourceLink final {
    SourceLink sourceLink;

    MapableSourceLink(SourceLink link) : sourceLink{std::move(link)} {}
  };
};

}  // namespace Acts

#include "Acts/AmbiguityResolution/GreedyAmbiguityResolution.ipp"

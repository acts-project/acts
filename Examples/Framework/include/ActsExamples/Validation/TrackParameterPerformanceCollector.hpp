// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/GeometryHierarchyMap.hpp"
#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/EventData/SimHit.hpp"
#include "ActsExamples/EventData/SimParticle.hpp"
#include "ActsExamples/EventData/Track.hpp"
#include "ActsExamples/EventData/TruthMatching.hpp"
#include "ActsExamples/Validation/ParametersOnSurface.hpp"
#include "ActsExamples/Validation/ResPlotTool.hpp"

#include <cstddef>
#include <memory>
#include <optional>
#include <vector>

namespace ActsExamples {

/// Collects residual/pull histograms for track-state parameters against the
/// truth hit information on the respective surface, without any file I/O.
///
/// The parameters of the selected track states are compared to the truth
/// obtained from the simulated hits associated to the state's measurement,
/// projected onto the state's reference surface. This works for track
/// parameters estimated from seeds, e.g. produced by `SeedsToTracks` which
/// stores the estimate as predicted parameters on the innermost track state,
/// as well as for fitted tracks.
///
/// @note The caller must ensure exclusive access (e.g. hold a mutex) when
///       calling fill(). This class applies no locking of its own.
class TrackParameterPerformanceCollector {
 public:
  struct Config {
    ResPlotTool::Config resPlotToolConfig = defaultResPlotToolConfig();
    /// The track-state parameters to compare to truth. If not set, the best
    /// available parameters are used (smoothed, filtered, or predicted).
    std::optional<TrackParameterType> parameterType;
    /// Optional geometry selection of track states. If non-empty, only track
    /// states within the given geometry hierarchy regions are used.
    std::vector<Acts::GeometryIdentifier> geometrySelection;
  };

  /// Default residual/pull plot configuration for bound parameters on a
  /// sensor surface, with binning suited for seed-estimated parameters.
  static ResPlotTool::Config defaultResPlotToolConfig();

  TrackParameterPerformanceCollector(
      Config cfg, std::unique_ptr<const Acts::Logger> logger);

  /// Fill histograms for one event.
  ///
  /// @note The caller must ensure exclusive access (e.g. hold a mutex).
  void fill(const Acts::GeometryContext& geoContext,
            const ConstTrackContainer& tracks,
            const SimParticleContainer& particles,
            const TrackParticleMatching& trackParticleMatching,
            const SimHitContainer& simHits,
            const MeasurementSimHitsMap& measurementSimHitsMap);

  /// Summary count statistics accumulated across all filled events.
  struct Stats {
    std::size_t nTotalTracks = 0;
    std::size_t nMatchedTracks = 0;
    std::size_t nFilledStates = 0;
  };

  /// Return accumulated event counts.
  const Stats& stats() const { return m_stats; }

  /// Emit summary statistics via the internal logger.
  void logSummary() const;

  /// Access to the underlying plot tool
  const ResPlotTool& resPlotTool() const { return m_resPlotTool; }

 private:
  const Acts::Logger& logger() const { return *m_logger; }

  Config m_cfg;
  std::unique_ptr<const Acts::Logger> m_logger;

  ResPlotTool m_resPlotTool;
  Acts::GeometryHierarchyMap<unsigned int> m_geometrySelection;

  Stats m_stats;
};

}  // namespace ActsExamples

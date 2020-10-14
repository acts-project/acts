// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/TrackFinding/CKFSourceLinkSelector.hpp"
#include "Acts/TrackFinding/CombinatorialKalmanFilter.hpp"
#include "ActsExamples/EventData/Track.hpp"
#include "ActsExamples/Framework/BareAlgorithm.hpp"
#include "ActsExamples/Plugins/BField/BFieldOptions.hpp"

#include <functional>
#include <vector>

namespace ActsExamples {

class TrackFindingAlgorithm final : public BareAlgorithm {
 public:
  /// Track finder function that takes input measurements, initial trackstate
  /// and track finder options and returns some track-finder-specific result.
  using TrackFinderOptions =
      Acts::CombinatorialKalmanFilterOptions<SimSourceLinkCalibrator,
                                             Acts::CKFSourceLinkSelector>;
  using TrackFinderResult =
      Acts::Result<Acts::CombinatorialKalmanFilterResult<SimSourceLink>>;
  using TrackFinderFunction = std::function<TrackFinderResult(
      const SimSourceLinkContainer&, const TrackParameters&,
      const TrackFinderOptions&)>;

  /// Create the track finder function implementation.
  ///
  /// The magnetic field is intentionally given by-value since the variant
  /// contains shared_ptr anyways.
  static TrackFinderFunction makeTrackFinderFunction(
      std::shared_ptr<const Acts::TrackingGeometry> trackingGeometry,
      Options::BFieldVariant magneticField);

  struct Config {
    /// Input source links collection.
    std::string inputSourceLinks;
    /// Input initial track parameter estimates for for each proto track.
    std::string inputInitialTrackParameters;
    /// Output find trajectories collection.
    std::string outputTrajectories;
    /// Type erased track finder function.
    TrackFinderFunction findTracks;
    /// CKF source link selector config
    Acts::CKFSourceLinkSelector::Config sourcelinkSelectorCfg;
  };

  /// Constructor of the track finding algorithm
  ///
  /// @param cfg is the config struct to configure the algorithm
  /// @param level is the logging level
  TrackFindingAlgorithm(Config cfg, Acts::Logging::Level lvl);

  /// Framework execute method of the track finding algorithm
  ///
  /// @param ctx is the algorithm context that holds event-wise information
  /// @return a process code to steer the algorithm flow
  ActsExamples::ProcessCode execute(
      const ActsExamples::AlgorithmContext& ctx) const final;

 private:
  Config m_cfg;
};

}  // namespace ActsExamples

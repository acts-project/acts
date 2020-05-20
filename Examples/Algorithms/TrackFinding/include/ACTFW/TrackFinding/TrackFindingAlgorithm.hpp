// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <functional>
//#include <memory>
#include <vector>

#include "ACTFW/EventData/SimSourceLink.hpp"
#include "ACTFW/EventData/Track.hpp"
#include "ACTFW/Framework/BareAlgorithm.hpp"
#include "ACTFW/Plugins/BField/BFieldOptions.hpp"
#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/TrackFinder/CKFSourceLinkSelector.hpp"
#include "Acts/TrackFinder/CombinatorialKalmanFilter.hpp"

namespace FW {

class TrackFindingAlgorithm final : public BareAlgorithm {
 public:
  using TrackFinderResult =
      Acts::Result<Acts::CombinatorialKalmanFilterResult<SimSourceLink>>;
  /// Track finding function that takes input measurements, initial trackstate
  /// and track finder options and returns some track-finding-specific result.
  using CKFOptions =
      Acts::CombinatorialKalmanFilterOptions<Acts::CKFSourceLinkSelector>;
  using TrackFinderFunction = std::function<TrackFinderResult(
      const SimSourceLinkContainer&, const TrackParameters&,
      const CKFOptions&)>;

  /// Create the track finder function implementation.
  ///
  /// The magnetic field is intentionally given by-value since the variant
  /// contains shared_ptr anyways.
  static TrackFinderFunction makeTrackFinderFunction(
      std::shared_ptr<const Acts::TrackingGeometry> trackingGeometry,
      Options::BFieldVariant magneticField, Acts::Logging::Level lvl);

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
  FW::ProcessCode execute(const FW::AlgorithmContext& ctx) const final override;

 private:
  Config m_cfg;
};

}  // namespace FW

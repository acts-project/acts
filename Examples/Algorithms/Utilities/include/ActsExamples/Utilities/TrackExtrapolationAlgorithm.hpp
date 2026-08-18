// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Utilities/Logger.hpp"
#include "Acts/Utilities/TrackHelpers.hpp"
#include "ActsExamples/EventData/Track.hpp"
#include "ActsExamples/Framework/DataHandle.hpp"
#include "ActsExamples/Framework/IAlgorithm.hpp"
#include "ActsExamples/Framework/ProcessCode.hpp"

#include <atomic>
#include <cstdint>
#include <memory>
#include <string>
#include <vector>

namespace Acts {
class MagneticFieldProvider;
class Surface;
class TrackingGeometry;
}  // namespace Acts

namespace ActsExamples {

/// Move the track parameters of a track container onto a common surface,
/// typically a perigee.
///
/// Only the track-level parameters are moved; the states are shared with the
/// input container and keep the parameters on their own surfaces. That is the
/// same layering a fitter produces.
///
/// Tracks whose extrapolation fails are dropped, so the output indices differ
/// from the input and any truth matching has to be redone downstream.
class TrackExtrapolationAlgorithm final : public IAlgorithm {
 public:
  struct Config {
    /// Input track collection.
    std::string inputTracks;
    /// Output track collection.
    std::string outputTracks;

    /// Surface to move the track parameters to.
    std::shared_ptr<const Acts::Surface> targetSurface;
    /// Tracking geometry to navigate.
    std::shared_ptr<const Acts::TrackingGeometry> trackingGeometry;
    /// Magnetic field to propagate in.
    std::shared_ptr<const Acts::MagneticFieldProvider> magneticField;

    /// Which track state to start the extrapolation from.
    Acts::TrackExtrapolationStrategy strategy =
        Acts::TrackExtrapolationStrategy::firstOrLast;

    /// The volume ids to constrain the extrapolation to.
    std::vector<std::uint32_t> constrainToVolumeIds;
    /// The volume ids to stop the extrapolation at.
    std::vector<std::uint32_t> endOfWorldVolumeIds;
  };

  /// Construct the algorithm.
  ///
  /// @param config the configuration
  /// @param logger the logger
  explicit TrackExtrapolationAlgorithm(
      Config config, std::unique_ptr<const Acts::Logger> logger = nullptr);

  /// Extrapolate the tracks of one event.
  ///
  /// @param ctx the algorithm context
  /// @return a process code
  ProcessCode execute(const AlgorithmContext& ctx) const override;

  /// Report how many tracks were dropped over the whole run.
  ///
  /// @return a process code
  ProcessCode finalize() override;

  /// Get readonly access to the config parameters
  /// @return the configuration
  const Config& config() const { return m_cfg; }

 private:
  Config m_cfg;

  mutable std::atomic<std::size_t> m_nTotalTracks{0};
  mutable std::atomic<std::size_t> m_nFailedTracks{0};

  ReadDataHandle<ConstTrackContainer> m_inputTracks{this, "InputTracks"};
  WriteDataHandle<ConstTrackContainer> m_outputTracks{this, "OutputTracks"};
};

}  // namespace ActsExamples

// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/EventData/Seed.hpp"
#include "ActsExamples/EventData/Track.hpp"
#include "ActsExamples/Framework/DataHandle.hpp"
#include "ActsExamples/Framework/IAlgorithm.hpp"
#include "ActsExamples/Framework/ProcessCode.hpp"

#include <memory>
#include <string>

namespace Acts {
class TrackingGeometry;
}

namespace ActsExamples {

/// Convert seeds into tracks with one track state per space-point source
/// link.
///
/// Optionally attaches estimated track parameters (parallel to the seeds) as
/// the track-level reference surface, parameters, and covariance, and sets
/// the track-state reference surfaces from the source-link geometry
/// identifiers.
class SeedsToTracks final : public IAlgorithm {
 public:
  struct Config {
    /// Input seeds.
    std::string inputSeeds = "seeds";
    /// Optional. Input track parameters, parallel to the input seeds, passed
    /// to the output tracks.
    std::string inputTrackParameters;
    /// Output tracks.
    std::string outputTracks = "tracks-from-seeds";
    /// Optional. If given, the reference surfaces of the track states are
    /// set from the source-link geometry identifiers.
    std::shared_ptr<const Acts::TrackingGeometry> trackingGeometry;
  };

  /// Construct the algorithm.
  ///
  /// @param cfg is the algorithm configuration
  /// @param logger is the logger
  explicit SeedsToTracks(Config cfg,
                         std::unique_ptr<const Acts::Logger> logger = nullptr);

  /// Run the algorithm.
  ///
  /// @param ctx is the algorithm context with event information
  /// @return a process code indication success or failure
  ProcessCode execute(const AlgorithmContext& ctx) const final;

  /// Const access to the config
  const Config& config() const { return m_cfg; }

 private:
  Config m_cfg;

  ReadDataHandle<SeedContainer> m_inputSeeds{this, "InputSeeds"};
  ReadDataHandle<TrackParametersContainer> m_inputTrackParameters{
      this, "InputTrackParameters"};
  WriteDataHandle<ConstTrackContainer> m_outputTracks{this, "OutputTracks"};
};

}  // namespace ActsExamples

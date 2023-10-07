// This file is part of the Acts project.
//
// Copyright (C) 2019-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "ActsExamples/EventData/ProtoTrack.hpp"
#include "ActsExamples/Framework/WriterT.hpp"

#include <memory>
#include <string>

namespace ActsExamples {

/// Write track finder performance measures.
///
/// Only considers the track finding itself, i.e. grouping of hits into tracks,
/// and computes relevant per-track and per-particles statistics.
class TrackFinderPerformanceWriter final : public WriterT<ProtoTrackContainer> {
 public:
  struct Config {
    /// Input reconstructed proto tracks collection.
    std::string inputProtoTracks;
    /// Input hit-particles map collection.
    std::string inputMeasurementParticlesMap;
    /// Input particles collection.
    std::string inputParticles;
    /// Output filename.
    std::string filePath = "performance_track_finder.root";
    /// Output file mode
    std::string fileMode = "RECREATE";
    /// Output tree name for the tracks
    std::string treeNameTracks = "track_finder_tracks";
    /// Output tree name for the particles
    std::string treeNameParticles = "track_finder_particles";
  };

  /// Constructor
  /// @param config the configuration
  /// @param level The log level
  TrackFinderPerformanceWriter(Config config, Acts::Logging::Level level);

  ~TrackFinderPerformanceWriter() override;

  ProcessCode finalize() override;

  /// Get readonly access to the config parameters
  const Config& config() const;

 private:
  ProcessCode writeT(const AlgorithmContext& ctx,
                     const ProtoTrackContainer& tracks) override;

  struct Impl;
  std::unique_ptr<Impl> m_impl;
};

}  // namespace ActsExamples

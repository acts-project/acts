// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <memory>
#include <string>

#include "ACTFW/EventData/ProtoTrack.hpp"
#include "ACTFW/Framework/WriterT.hpp"

namespace FW {

/// Write track finder performance measures.
///
/// Only considers the track finding itself, i.e. grouping of hits into tracks,
/// and computes relevant per-track and per-particles statistics.
class TrackFinderPerformanceWriter final : public WriterT<ProtoTrackContainer> {
 public:
  struct Config {
    /// True set of input particles.
    std::string inputParticles;
    /// True hit-particles mapping.
    std::string inputHitParticlesMap;
    /// Reconstructed input proto tracks.
    std::string inputProtoTracks;
    /// Output directory.
    std::string outputDir;
    /// Output filename
    std::string outputFilename = "performance_track_finder.root";
  };

  TrackFinderPerformanceWriter(Config cfg, Acts::Logging::Level lvl);
  ~TrackFinderPerformanceWriter();

  ProcessCode endRun() final override;

 private:
  ProcessCode writeT(const AlgorithmContext& ctx,
                     const ProtoTrackContainer& tracks) final override;

  struct Impl;
  std::unique_ptr<Impl> m_impl;
};

}  // namespace FW

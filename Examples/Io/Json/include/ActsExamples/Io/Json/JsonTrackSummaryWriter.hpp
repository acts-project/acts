// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/EventData/Jets.hpp"
#include "ActsExamples/EventData/SimParticle.hpp"
#include "ActsExamples/EventData/Track.hpp"
#include "ActsExamples/EventData/TruthMatching.hpp"
#include "ActsExamples/Framework/DataHandle.hpp"
#include "ActsExamples/Framework/ProcessCode.hpp"
#include "ActsExamples/Framework/WriterT.hpp"

#include <mutex>
#include <string>
#include <vector>

#include <nlohmann/json.hpp>

namespace ActsExamples {

/// Write out reconstructed track summary information as JSON.
///
/// Mirrors RootTrackSummaryWriter exactly: writes the same per-track fields
/// (track quality counters, fitted parameters, errors, residuals, pulls, truth
/// matching) without depending on ROOT.  The output is a JSON array of event
/// objects; each event contains an array of track objects.
///
/// Safe to use from multiple writer threads - uses a std::mutex lock.
class JsonTrackSummaryWriter final : public WriterT<ConstTrackContainer> {
 public:
  struct Config {
    /// Input (fitted) tracks collection.
    std::string inputTracks;
    /// Input particles collection (optional).
    std::string inputParticles;
    /// Input track-particle matching (optional).
    std::string inputTrackParticleMatching;
    /// Output filename.
    std::string filePath = "tracksummary.json";
    /// Write full 6x6 covariance matrix entries.
    bool writeCovMat = false;
    /// Write GSF material statistics (kFwdMaxMaterialXOverX0 /
    /// kFwdSumMaterialXOverX0).
    bool writeGsfSpecific = false;
    /// Write GX2F update count.
    bool writeGx2fSpecific = false;
  };

  /// Constructor
  JsonTrackSummaryWriter(const Config& config, Acts::Logging::Level level);
  ~JsonTrackSummaryWriter() override;

  /// Write accumulated events to the JSON output file.
  ProcessCode finalize() override;

  const Config& config() const { return m_cfg; }

 protected:
  ProcessCode writeT(const AlgorithmContext& ctx,
                     const ConstTrackContainer& tracks) override;

 private:
  Config m_cfg;

  ReadDataHandle<SimParticleContainer> m_inputParticles{this, "InputParticles"};
  ReadDataHandle<TrackParticleMatching> m_inputTrackParticleMatching{
      this, "InputTrackParticleMatching"};

  std::mutex m_writeMutex;
  /// One JSON object per event, accumulated across all calls to writeT.
  std::vector<nlohmann::json> m_events;
};

}  // namespace ActsExamples

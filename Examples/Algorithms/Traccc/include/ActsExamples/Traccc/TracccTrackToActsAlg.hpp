// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/EventData/Track.hpp"
#include "ActsExamples/Framework/DataHandle.hpp"
#include "ActsExamples/Framework/IAlgorithm.hpp"
#include "ActsExamples/Traccc/TracccDetector.hpp"
#include "ActsExamples/Traccc/TracccEventData.hpp"

#include <memory>
#include <string>


namespace ActsExamples {

/// Convert (host) traccc tracks into ACTS tracks
///
/// The track states are linked to the ACTS measurements through the
/// identifier written by the @c ActsMeasToTracccAlg, so the resulting tracks
/// can be used with the standard ACTS truth matching and performance writers.
class TracccTrackToActsAlg final : public IAlgorithm {
 public:
  struct Config {
    /// The traccc detector, providing the backend and surface mapping
    std::shared_ptr<const TracccDetector> detector;
    /// Input host traccc tracks
    std::string inputTracks = "traccc_host_tracks";
    /// Input host copy of the (sorted) traccc measurements the tracks refer to
    std::string inputMeasurements = "traccc_host_measurements_sorted";
    /// Output ACTS tracks
    std::string outputTracks = "traccc_tracks";
  };

  explicit TracccTrackToActsAlg(
      const Config& cfg, std::unique_ptr<const Acts::Logger> logger = nullptr);

  ProcessCode execute(const AlgorithmContext& ctx) const override;

  const Config& config() const { return m_cfg; }

 private:
  Config m_cfg;

  ReadDataHandle<TracccHostTracks> m_inputTracks{this, "InputTracks"};
  ReadDataHandle<TracccHostMeasurements> m_inputMeasurements{
      this, "InputMeasurements"};
  WriteDataHandle<ConstTrackContainer> m_outputTracks{this, "OutputTracks"};
};

}  // namespace ActsExamples

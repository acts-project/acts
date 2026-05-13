// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "ActsExamples/EventData/Measurement.hpp"
#include "ActsExamples/EventData/ProtoTrack.hpp"
#include "ActsExamples/EventData/Track.hpp"
#include "ActsExamples/Framework/DataHandle.hpp"
#include "ActsExamples/Framework/IAlgorithm.hpp"

namespace ActsExamples {

class TracccTracksToActsTracks final : public IAlgorithm {
 public:
  struct Config {
    /// Input proto tracks.
    std::string inputProtoTracks;
    /// Optional. Input track parameters passed to the output tracks.
    std::string inputTrackParameters;
    /// Input measurements.
    std::string inputMeasurements;
    /// Output tracks.
    std::string outputTracks = "tracks-from-tracccTracks";
  };

  /// Construct the algorithm.
  ///
  /// @param cfg is the algorithm configuration
  /// @param lvl is the logging level
  explicit TracccTracksToActsTracks(
      Config cfg, std::unique_ptr<const Acts::Logger> logger = nullptr);
  /// Run the algorithm.
  ///
  /// @param ctx is the algorithm context with event information
  /// @return a process code indication success or failure
  ProcessCode execute(const AlgorithmContext& ctx) const final;

  /// Const access to the config
  const Config& config() const { return m_cfg; }

 private:
  Config m_cfg;

  WriteDataHandle<ConstTrackContainer> m_outputTracks{this, "OutputActsTracks"};
  ReadDataHandle<MeasurementSubset> m_inputMeasurements{this,
                                                        "InputTracccMeasurements"};
  ReadDataHandle<TrackParametersContainer> m_inputTrackParameters{
      this, "InputTracccTrackParameters"};
  ReadDataHandle<ProtoTrackContainer> m_inputProtoTracks{this,
                                                         "InputTracccTracks"};
};

}  // namespace ActsExamples
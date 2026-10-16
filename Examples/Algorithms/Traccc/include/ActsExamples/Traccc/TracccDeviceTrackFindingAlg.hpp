// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/Framework/DataHandle.hpp"
#include "ActsExamples/Framework/IAlgorithm.hpp"
#include "ActsExamples/Traccc/TracccBackend.hpp"
#include "ActsExamples/Traccc/TracccDetector.hpp"
#include "ActsExamples/Traccc/TracccEventData.hpp"
#include "ActsExamples/Traccc/TracccMagneticField.hpp"

#include <memory>
#include <optional>
#include <string>

#include <traccc/finding/finding_config.hpp>

namespace ActsExamples {

/// Run the traccc combinatorial Kalman filter on the device
class TracccDeviceTrackFindingAlg final : public IAlgorithm {
 public:
  struct Config {
    /// The traccc detector, providing the backend and geometry
    std::shared_ptr<const TracccDetector> detector;
    /// The traccc magnetic field
    std::shared_ptr<const TracccMagneticField> magneticField;
    /// Input sorted device traccc measurements
    std::string inputMeasurements = "traccc_device_measurements";
    /// Input device traccc track parameters
    std::string inputTrackParameters = "traccc_device_track_parameters";
    /// Output device traccc tracks
    std::string outputTracks = "traccc_device_tracks";
    /// Track finding configuration
    traccc::finding_config finding;
  };

  explicit TracccDeviceTrackFindingAlg(
      const Config& cfg, std::unique_ptr<const Acts::Logger> logger = nullptr);

  ProcessCode execute(const AlgorithmContext& ctx) const override;

  const Config& config() const { return m_cfg; }

 private:
  Config m_cfg;
  std::optional<TracccDeviceAlgorithm<
      traccc::device::combinatorial_kalman_filter_algorithm>>
      m_algorithm;

  ReadDataHandle<TracccDeviceMeasurements> m_inputMeasurements{
      this, "InputMeasurements"};
  ReadDataHandle<TracccDeviceTrackParameters> m_inputTrackParameters{
      this, "InputTrackParameters"};
  WriteDataHandle<TracccDeviceTracks> m_outputTracks{this, "OutputTracks"};
};

}  // namespace ActsExamples

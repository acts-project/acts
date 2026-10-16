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

#include <traccc/seeding/detail/track_params_estimation_config.hpp>

namespace ActsExamples {

/// Estimate the track parameters of traccc seeds on the device
class TracccDeviceTrackParamEstimationAlg final : public IAlgorithm {
 public:
  struct Config {
    /// The traccc detector, providing the backend
    std::shared_ptr<const TracccDetector> detector;
    /// The traccc magnetic field
    std::shared_ptr<const TracccMagneticField> magneticField;
    /// Input device traccc measurements
    std::string inputMeasurements = "traccc_device_measurements";
    /// Input device traccc spacepoints
    std::string inputSpacepoints = "traccc_device_spacepoints";
    /// Input device traccc seeds
    std::string inputSeeds = "traccc_device_seeds";
    /// Output device traccc track parameters
    std::string outputTrackParameters = "traccc_device_track_parameters";
    /// Parameter estimation configuration
    traccc::track_params_estimation_config estimation;
  };

  explicit TracccDeviceTrackParamEstimationAlg(
      const Config& cfg, std::unique_ptr<const Acts::Logger> logger = nullptr);

  ProcessCode execute(const AlgorithmContext& ctx) const override;

  const Config& config() const { return m_cfg; }

 private:
  Config m_cfg;
  std::optional<TracccDeviceAlgorithm<
      traccc::device::seed_parameter_estimation_algorithm>>
      m_algorithm;

  ReadDataHandle<TracccDeviceMeasurements> m_inputMeasurements{
      this, "InputMeasurements"};
  ReadDataHandle<TracccDeviceSpacepoints> m_inputSpacepoints{
      this, "InputSpacepoints"};
  ReadDataHandle<TracccDeviceSeeds> m_inputSeeds{this, "InputSeeds"};
  WriteDataHandle<TracccDeviceTrackParameters> m_outputTrackParameters{
      this, "OutputTrackParameters"};
};

}  // namespace ActsExamples

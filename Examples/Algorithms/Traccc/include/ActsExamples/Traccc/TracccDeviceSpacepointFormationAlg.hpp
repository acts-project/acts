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

#include <memory>
#include <optional>
#include <string>

namespace ActsExamples {

/// Form spacepoints from the 2D (pixel) measurements on the device
class TracccDeviceSpacepointFormationAlg final : public IAlgorithm {
 public:
  struct Config {
    /// The traccc detector, providing the backend and the device geometry
    std::shared_ptr<const TracccDetector> detector;
    /// Input device traccc measurements
    std::string inputMeasurements = "traccc_device_measurements";
    /// Output device traccc spacepoints
    std::string outputSpacepoints = "traccc_device_spacepoints";
  };

  explicit TracccDeviceSpacepointFormationAlg(
      const Config& cfg, std::unique_ptr<const Acts::Logger> logger = nullptr);

  ProcessCode execute(const AlgorithmContext& ctx) const override;

  const Config& config() const { return m_cfg; }

 private:
  Config m_cfg;
  std::optional<TracccDeviceAlgorithm<
      traccc::device::silicon_pixel_spacepoint_formation_algorithm>>
      m_algorithm;

  ReadDataHandle<TracccDeviceMeasurements> m_inputMeasurements{
      this, "InputMeasurements"};
  WriteDataHandle<TracccDeviceSpacepoints> m_outputSpacepoints{
      this, "OutputSpacepoints"};
};

}  // namespace ActsExamples

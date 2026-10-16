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

/// Sort device traccc measurements by surface, as required by the traccc
/// track finding
class TracccDeviceMeasurementSortingAlg final : public IAlgorithm {
 public:
  struct Config {
    /// The traccc detector, providing the backend
    std::shared_ptr<const TracccDetector> detector;
    /// Input device traccc measurements
    std::string inputMeasurements = "traccc_device_measurements_unsorted";
    /// Output sorted device traccc measurements
    std::string outputMeasurements = "traccc_device_measurements";
  };

  explicit TracccDeviceMeasurementSortingAlg(
      const Config& cfg, std::unique_ptr<const Acts::Logger> logger = nullptr);

  ProcessCode execute(const AlgorithmContext& ctx) const override;

  const Config& config() const { return m_cfg; }

 private:
  Config m_cfg;
  std::optional<
      TracccDeviceAlgorithm<TracccBackend::MeasurementSortingAlgorithm>>
      m_sorting;

  ReadDataHandle<TracccDeviceMeasurements> m_inputMeasurements{
      this, "InputMeasurements"};
  WriteDataHandle<TracccDeviceMeasurements> m_outputMeasurements{
      this, "OutputMeasurements"};
};

}  // namespace ActsExamples

// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/EventData/Measurement.hpp"
#include "ActsExamples/Framework/DataHandle.hpp"
#include "ActsExamples/Framework/IAlgorithm.hpp"
#include "ActsExamples/Traccc/TracccDetector.hpp"
#include "ActsExamples/Traccc/TracccEventData.hpp"

#include <memory>
#include <string>

#include <vecmem/memory/host_memory_resource.hpp>

namespace ActsExamples {

/// Convert ACTS measurements into a host traccc measurement collection
///
/// The ACTS measurement index is stored as the traccc measurement identifier,
/// so that tracks found by traccc can be linked back to the ACTS measurements.
/// Only the local position components are converted. 2D measurements become
/// spacepoints in the traccc spacepoint formation.
class ActsMeasToTracccAlg final : public IAlgorithm {
 public:
  struct Config {
    /// The traccc detector, providing the ACTS to detray surface mapping
    std::shared_ptr<const TracccDetector> detector;
    /// Input ACTS measurements
    std::string inputMeasurements = "measurements";
    /// Output host traccc measurements
    std::string outputMeasurements = "traccc_host_measurements";
  };

  explicit ActsMeasToTracccAlg(
      const Config& cfg, std::unique_ptr<const Acts::Logger> logger = nullptr);

  ProcessCode execute(const AlgorithmContext& ctx) const override;

  const Config& config() const { return m_cfg; }

 private:
  Config m_cfg;
  mutable vecmem::host_memory_resource m_mr;

  ReadDataHandle<MeasurementContainer> m_inputMeasurements{this,
                                                           "InputMeasurements"};
  WriteDataHandle<TracccHostMeasurements> m_outputMeasurements{
      this, "OutputMeasurements"};
};

}  // namespace ActsExamples

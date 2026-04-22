// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "ActsExamples/EventData/Measurement.hpp"
#include "ActsExamples/Framework/DataHandle.hpp"
#include "ActsExamples/Framework/IAlgorithm.hpp"
#include "ActsExamples/Framework/ProcessCode.hpp"

#include <string>

namespace ActsExamples {

/// Produces a MeasurementContainer that excludes measurements already used
/// by a previous tracking pass.
///
/// This is the standalone-Acts equivalent of the Athena cluster-preparation
/// and spacepoint-preparation algorithms that filter clusters/spacepoints
/// using a PrepRawDataAssociation map before the next tracking pass begins.
///
/// The output container is a compact copy: measurements are re-indexed
/// 0..M-1 where M is the number of unused measurements.  The accompanying
/// outputIndexRemapping vector maps each new index back to the corresponding
/// index in the original container, allowing MeasurementMapAlgorithm in the
/// next pass to convert track-state source-link indices back to
/// original-container space when accumulating the map.
class MeasurementFilterAlgorithm final : public IAlgorithm {
 public:
  struct Config {
    /// Input MeasurementContainer (original, from digitization).
    std::string inputMeasurements;
    /// MeasurementMap of already-used measurement indices (original-container
    /// space) produced by MeasurementMapAlgorithm from the previous pass.
    std::string inputMeasurementMap;
    /// Output filtered MeasurementContainer (re-indexed 0..M-1).
    std::string outputMeasurements;
    /// Output index remapping: element i = original index of output
    /// measurement i.  Consumed by MeasurementMapAlgorithm.
    std::string outputIndexRemapping;
  };

  explicit MeasurementFilterAlgorithm(
      Config cfg, std::unique_ptr<const Acts::Logger> logger = nullptr);

  ProcessCode execute(const AlgorithmContext& ctx) const final;

  const Config& config() const { return m_cfg; }

 private:
  Config m_cfg;

  ReadDataHandle<MeasurementContainer> m_inputMeasurements{this,
                                                           "InputMeasurements"};
  ReadDataHandle<UsedMeasurementMap> m_inputMeasurementMap{
      this, "InputMeasurementMap"};
  WriteDataHandle<MeasurementContainer> m_outputMeasurements{
      this, "OutputMeasurements"};
  WriteDataHandle<MeasurementIndexRemapping> m_outputIndexRemapping{
      this, "OutputIndexRemapping"};
};

}  // namespace ActsExamples

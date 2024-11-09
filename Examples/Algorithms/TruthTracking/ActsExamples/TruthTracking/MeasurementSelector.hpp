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
#include "ActsExamples/EventData/SimParticle.hpp"
#include "ActsExamples/Framework/DataHandle.hpp"
#include "ActsExamples/Framework/IAlgorithm.hpp"

#include <limits>
#include <string>

namespace ActsExamples {

/// Select measurements by applying some selection cuts.
class MeasurementSelector final : public IAlgorithm {
 public:
  struct Config {
    /// Input measurement collection.
    std::string inputMeasurements;
    /// Input measurement particles map.
    std::string inputMeasurementParticlesMap;
    /// Optional input particle collection.
    std::string inputParticlesSelected;
    /// Output measurement collection
    std::string outputMeasurements;
    /// Output measurement particles map
    std::string outputMeasurementParticlesMap;

    /// Min primary vertex ID cut
    std::uint64_t minPrimaryVertexId = 0;
    /// Max primary vertex ID cut
    std::uint64_t maxPrimaryVertexId =
        std::numeric_limits<std::uint64_t>::max();
  };

  MeasurementSelector(const Config& config, Acts::Logging::Level level);

  ProcessCode execute(const AlgorithmContext& ctx) const final;

  /// Get readonly access to the config parameters
  const Config& config() const { return m_cfg; }

 private:
  Config m_cfg;

  ReadDataHandle<MeasurementContainer> m_inputMeasurements{this,
                                                           "InputMeasurements"};
  ReadDataHandle<IndexMultimap<SimBarcode>> m_inputMeasurementParticlesMap{
      this, "MeasurementParticlesMap"};
  ReadDataHandle<SimParticleContainer> m_inputParticlesSelected{
      this, "InputParticlesSelected"};
  WriteDataHandle<MeasurementContainer> m_outputMeasurements{
      this, "OutputMeasurements"};
  WriteDataHandle<IndexMultimap<SimBarcode>> m_outputMeasurementParticlesMap{
      this, "OutputMeasurementParticlesMap"};
};

}  // namespace ActsExamples

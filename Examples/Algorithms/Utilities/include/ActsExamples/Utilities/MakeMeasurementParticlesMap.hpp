// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "ActsExamples/EventData/Index.hpp"
#include "ActsExamples/EventData/Measurement.hpp"
#include "ActsExamples/EventData/SimHit.hpp"
#include "ActsExamples/Framework/DataHandle.hpp"
#include "ActsExamples/Framework/IAlgorithm.hpp"

namespace ActsExamples {

class MakeMeasurementParticlesMap final : public IAlgorithm {
 public:
  struct Config {
    std::string inputSimHits;
    std::string inputMeasurementSimhitMap;
    std::string outputMeasurementParticlesMap;
  };

  /// Construct the algorithm.
  ///
  /// @param cfg is the algorithm configuration
  /// @param lvl is the logging level
  MakeMeasurementParticlesMap(Config cfg, Acts::Logging::Level lvl);

  /// Run the algorithm.
  ///
  /// @param ctx is the algorithm context with event information
  /// @return a process code indication success or failure
  ProcessCode execute(const AlgorithmContext& ctx) const final;

  /// Const access to the config
  const Config& config() const { return m_cfg; }

 private:
  Config m_cfg;

  WriteDataHandle<IndexMultimap<ActsFatras::Barcode>> m_outputParticleMap{
      this, "OutputMeasurementParticlesMap"};
  ReadDataHandle<IndexMultimap<Index>> m_inputHitMap{
      this, "InputMeasurementSimhitMap"};
  ReadDataHandle<SimHitContainer> m_inputHits{this, "InputHits"};
};

}  // namespace ActsExamples

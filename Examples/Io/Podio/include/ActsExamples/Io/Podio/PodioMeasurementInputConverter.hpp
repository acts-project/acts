// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "ActsExamples/EventData/Measurement.hpp"
#include "ActsExamples/EventData/SimHit.hpp"
#include "ActsExamples/Framework/DataHandle.hpp"
#include "ActsExamples/Io/Podio/PodioInputConverter.hpp"
#include "ActsPlugins/EDM4hep/EDM4hepUtil.hpp"

namespace ActsExamples {

class PodioMeasurementInputConverter : public PodioInputConverter {
 public:
  struct Config {
    std::string inputFrame;
    /// Input measurement collection name in podio
    std::string inputMeasurements = "ActsMeasurements";
    /// Output measurement collection
    std::string outputMeasurements;
    /// Output collection to map measured hits to contributing particles.
    std::string outputMeasurementParticlesMap;
    /// Output collection to map measured hits to simulated hits.
    std::string outputMeasurementSimHitsMap;
    /// Output collection to map particles to measurements.
    std::string outputParticleMeasurementsMap;
    /// Output collection to map simulated hits to measurements.
    std::string outputSimHitMeasurementsMap;

    /// Sim hit collection used to associated read measurements
    std::string inputSimHits;
    /// Sim hit association used to connect measurements to sim hits
    std::string inputSimHitAssociation;
  };

  /// constructor
  /// @param config is the configuration object
  /// @param level is the output logging level
  explicit PodioMeasurementInputConverter(
      const Config& config, Acts::Logging::Level level = Acts::Logging::INFO);

  /// Readonly access to the config
  const Config& config() const { return m_cfg; }

  ProcessCode convert(const AlgorithmContext& ctx,
                      const podio::Frame& frame) const final;

 private:
  Config m_cfg;

  WriteDataHandle<MeasurementContainer> m_outputMeasurements{
      this, "OutputMeasurements"};

  WriteDataHandle<IndexMultimap<SimBarcode>> m_outputMeasurementParticlesMap{
      this, "OutputMeasurementParticlesMap"};

  WriteDataHandle<IndexMultimap<Index>> m_outputMeasurementSimHitsMap{
      this, "OutputMeasurementSimHitsMap"};

  WriteDataHandle<InverseMultimap<SimBarcode>> m_outputParticleMeasurementsMap{
      this, "OutputParticleMeasurementsMap"};

  WriteDataHandle<InverseMultimap<Index>> m_outputSimHitMeasurementsMap{
      this, "OutputSimHitMeasurementsMap"};

  // Temporary!
  ReadDataHandle<SimHitContainer> m_inputSimHits{this, "InputSimHits"};
  ReadDataHandle<ActsPlugins::EDM4hepUtil::SimHitAssociation>
      m_inputSimHitAssociation{this, "InputSimHitAssociation"};
};

}  // namespace ActsExamples

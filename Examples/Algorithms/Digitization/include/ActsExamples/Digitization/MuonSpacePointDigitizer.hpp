// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "Acts/Utilities/RangeXD.hpp"
#include "ActsExamples/EventData/Measurement.hpp"
#include "ActsExamples/EventData/MuonSpacePoint.hpp"
#include "ActsExamples/EventData/MuonSpacePointCalibrator.hpp"
#include "ActsExamples/EventData/SimHit.hpp"
#include "ActsExamples/EventData/SimParticle.hpp"
#include "ActsExamples/Framework/DataHandle.hpp"
#include "ActsExamples/Framework/IAlgorithm.hpp"
#include "ActsExamples/Framework/ProcessCode.hpp"
#include "ActsExamples/Framework/RandomNumbers.hpp"

namespace ActsExamples {

/// Algorithm that turns simulated hits into measurements by truth smearing.
class MuonSpacePointDigitizer final : public IAlgorithm {
 public:
  struct Config {
    /// @brief Name of the input simulated hits collection
    std::string inputSimHits{"simhits"};
    /// @brief Name of the input simulated particles collection
    std::string inputParticles{"particles_simulated"};
    /// @brief Name of the output spacepoints collection
    std::string outputSpacePoints{"MuonSpacePoints"};
    /// @brief Output measurements collection.
    std::string outputMeasurements = "measurements";
    /// @brief Output collection to map measured hits to contributing particles.
    std::string outputMeasurementParticlesMap = "measurement_particles_map";
    /// @brief Output collection to map measured hits to simulated hits.
    std::string outputMeasurementSimHitsMap = "measurement_simhits_map";
    /// @brief Output collection to map particles to measurements.
    std::string outputParticleMeasurementsMap = "particle_measurements_map";
    /// @brief Output collection to map particles to simulated hits.
    std::string outputSimHitMeasurementsMap = "simhit_measurements_map";
    /// @brief Random number generator service
    std::shared_ptr<const RandomNumbers> randomNumbers{};
    /// @brief Pointer to the muon calibrator to fetch the smearing constants
    std::shared_ptr<const MuonSpacePointCalibrator> calibrator{};
    /// @brief Pointer to the tracking geometry to fetch the surfaces
    std::shared_ptr<const Acts::TrackingGeometry> trackingGeometry{};
    /// @brief Add strip times
    bool digitizeTime{false};
    /// @brief Visualize the digitization
    bool dumpVisualization{true};
    /// @brief Visualization function (optional, e.g., for ROOT-based visualization)
    /// Takes: outputPath, gctx, bucket, simHits, simParticles,
    /// trackingGeometry, logger
    std::function<void(const std::string&, const Acts::GeometryContext&,
                       const MuonSpacePointBucket&, const SimHitContainer&,
                       const SimParticleContainer&,
                       const Acts::TrackingGeometry&, const Acts::Logger&)>
        visualizationFunction{};
    /// @brief Applied dead time between two consecutive straw hits
    double strawDeadTime{1. * Acts::UnitConstants::ms};
    /// @brief Applied dead time between two consecutive rpc hits
    double rpcDeadTime{50. * Acts::UnitConstants::ns};
  };
  /// @brief Constructor
  MuonSpacePointDigitizer(const Config& cfg, Acts::Logging::Level lvl);

  /// @brief Initialize the digitizer
  ProcessCode initialize() override;
  /// @brief Execute the digitization
  ProcessCode execute(const AlgorithmContext& ctx) const override;
  /// @brief  Const access to the config
  const Config& config() const { return m_cfg; }

 private:
  /// @brief Returns the reference to the configured calibrator
  const MuonSpacePointCalibrator& calibrator() const;
  /// @brief Returns the reference to the passed tracking geometry
  const Acts::TrackingGeometry& trackingGeometry() const;
  /// @brief  Returns the transformation from the local hit frame into the
  ///         chamber's surface frame
  /// @param gctx Geometry context to access the local -> global transform of the surface
  /// @param hitId Geometry identifier of the hit of interest
  Acts::Transform3 toSpacePointFrame(
      const Acts::GeometryContext& gctx,
      const Acts::GeometryIdentifier& hitId) const;
  /// @brief Configuration of the digitizer
  Config m_cfg;
  /// @brief Data handle for the input simulated hits
  ReadDataHandle<SimHitContainer> m_inputSimHits{this, "InputSimHits"};
  /// @brief Data handle for the input simulated particles
  ReadDataHandle<SimParticleContainer> m_inputParticles{this, "G4Particles"};
  /// @brief Data handle for the output spacepoints
  WriteDataHandle<MuonSpacePointContainer> m_outputSpacePoints{this,
                                                               "SpacePoints"};
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
};
}  // namespace ActsExamples

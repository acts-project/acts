// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Utilities/Histogram.hpp"
#include "ActsExamples/EventData/Measurement.hpp"
#include "ActsExamples/EventData/SpacePoint.hpp"
#include "ActsExamples/EventData/TruthMatching.hpp"
#include "ActsExamples/Framework/DataHandle.hpp"
#include "ActsExamples/Framework/ProcessCode.hpp"
#include "ActsExamples/Framework/WriterT.hpp"
#include "ActsExamples/Utilities/StripModulePairing.hpp"

#include <string>

class TFile;
class TTree;

namespace ActsExamples {

class RootSpacePointPerformanceWriter final
    : public WriterT<SpacePointContainer> {
 public:
  struct Config {
    /// Which space point collection to write.
    std::string inputSpacePoints;
    /// Input particle collection for truth matching.
    std::string inputParticles;
    /// Which measurement collection to write.
    std::string inputMeasurements;
    /// Which simulated (truth) hits collection to use.
    std::string inputSimHits;
    /// Input collection to map reconstructed measurements to simulated hits.
    std::string inputMeasurementSimHitsMap;
    /// Input collection to map reconstructed measurements to particles.
    std::string inputMeasurementParticlesMap;

    /// Tracking geometry for transformation lookup.
    std::shared_ptr<const Acts::TrackingGeometry> trackingGeometry;
    /// Geometry selection for strip modules
    std::vector<Acts::GeometryIdentifier> stripGeometrySelection;

    /// path of the output file
    std::string filePath = "";
    /// file access mode
    std::string fileMode = "RECREATE";
    /// The tree name
    std::string treeName = "spacePoints";

    Acts::Experimental::AxisVariant zAxis =
        Acts::Experimental::BoostRegularAxis{100, -3000, 3000, "z [mm]"};
    Acts::Experimental::AxisVariant rAxis =
        Acts::Experimental::BoostRegularAxis{100, 0, 1200, "r [mm]"};
    Acts::Experimental::AxisVariant etaAxis =
        Acts::Experimental::BoostRegularAxis{40, -3, 3, "#eta"};
    Acts::Experimental::AxisVariant phiAxis =
        Acts::Experimental::BoostRegularAxis{100, -std::numbers::pi,
                                             std::numbers::pi, "#phi"};
  };

  /// Constructor with
  /// @param cfg configuration struct
  /// @param output logging level
  RootSpacePointPerformanceWriter(const Config& config,
                                  Acts::Logging::Level level);

  /// Virtual destructor
  ~RootSpacePointPerformanceWriter() override;

  /// End-of-run hook
  ProcessCode finalize() override;

  /// Get const access to the config
  const Config& config() const { return m_cfg; }

 protected:
  /// This implementation holds the actual writing method
  /// and is called by the WriterT<>::write interface
  ///
  /// @param ctx The Algorithm context with per event information
  /// @param spacePoints is the data to be written out
  ProcessCode writeT(const AlgorithmContext& ctx,
                     const SpacePointContainer& spacePoints) override;

 private:
  Config m_cfg;

  StripModulePairMap m_stripModulePairMap;

  ReadDataHandle<SimParticleContainer> m_inputParticles{this, "InputParticles"};
  ReadDataHandle<MeasurementContainer> m_inputMeasurements{this,
                                                           "InputMeasurements"};
  ReadDataHandle<SimHitContainer> m_inputSimHits{this, "InputSimHits"};

  ReadDataHandle<MeasurementSimHitsMap> m_inputMeasurementSimHitsMap{
      this, "InputMeasurementSimHitsMap"};
  ReadDataHandle<MeasurementParticlesMap> m_inputMeasurementParticlesMap{
      this, "InputMeasurementParticlesMap"};

  /// protect multi-threaded writes
  std::mutex m_writeMutex;
  /// the output file
  TFile* m_outputFile = nullptr;

  std::optional<Acts::Experimental::Efficiency1> m_fakeVsZ;
  std::optional<Acts::Experimental::Efficiency1> m_fakeVsR;
  std::optional<Acts::Experimental::Efficiency1> m_fakeVsEta;
  std::optional<Acts::Experimental::Efficiency1> m_fakeVsPhi;

  std::optional<Acts::Experimental::Efficiency1> m_effVsZ;
  std::optional<Acts::Experimental::Efficiency1> m_effVsR;
  std::optional<Acts::Experimental::Efficiency1> m_effVsEta;
  std::optional<Acts::Experimental::Efficiency1> m_effVsPhi;
};

}  // namespace ActsExamples

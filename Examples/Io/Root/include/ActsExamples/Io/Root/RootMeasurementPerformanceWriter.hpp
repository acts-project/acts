// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Utilities/Histogram.hpp"
#include "ActsExamples/EventData/Index.hpp"
#include "ActsExamples/EventData/Measurement.hpp"
#include "ActsExamples/EventData/SimHit.hpp"
#include "ActsExamples/Framework/DataHandle.hpp"
#include "ActsExamples/Framework/ProcessCode.hpp"
#include "ActsExamples/Framework/WriterT.hpp"

#include <string>

class TFile;
class TTree;

namespace ActsExamples {

class RootMeasurementPerformanceWriter final
    : public WriterT<MeasurementContainer> {
 public:
  enum class MeasurementClassification {
    Unknown = 0,
    Matched,
    Merged,
    Fake,
  };

  struct Config {
    /// Which measurement collection to write.
    std::string inputMeasurements;
    /// Which simulated (truth) hits collection to use.
    std::string inputSimHits;
    /// Input collection to map reconstructed measurements to simulated hits.
    std::string inputMeasurementSimHitsMap;
    /// Input collection to map reconstructed measurements to particles.
    std::string inputMeasurementParticlesMap;
    /// Input collection to map simulated hits to reconstructed measurements.
    std::string inputSimHitMeasurementsMap;

    /// path of the output file
    std::string filePath = "";
    /// file access mode
    std::string fileMode = "RECREATE";
    /// The tree name
    std::string treeName = "measurements";

    /// Matching ratio threshold for measurement to particle matching
    double matchingRatio = 0.5;

    Acts::Experimental::AxisVariant countAxis =
        Acts::Experimental::BoostRegularAxis{10, 0, 10, "Count"};
    Acts::Experimental::AxisVariant purityAxis =
        Acts::Experimental::BoostRegularAxis{10, 0, 1.1, "Purity"};
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
  RootMeasurementPerformanceWriter(const Config& config,
                                   Acts::Logging::Level level);

  /// Virtual destructor
  ~RootMeasurementPerformanceWriter() override;

  /// End-of-run hook
  ProcessCode finalize() override;

  /// Get const access to the config
  const Config& config() const { return m_cfg; }

 protected:
  /// This implementation holds the actual writing method
  /// and is called by the WriterT<>::write interface
  ///
  /// @param ctx The Algorithm context with per event information
  /// @param measurements is the data to be written out
  ProcessCode writeT(const AlgorithmContext& ctx,
                     const MeasurementContainer& measurements) override;

 private:
  Config m_cfg;

  ReadDataHandle<SimHitContainer> m_inputSimHits{this, "InputSimHits"};

  ReadDataHandle<IndexMultimap<Index>> m_inputMeasurementSimHitsMap{
      this, "InputMeasurementSimHitsMap"};
  ReadDataHandle<IndexMultimap<SimBarcode>> m_inputMeasurementParticlesMap{
      this, "InputMeasurementParticlesMap"};
  ReadDataHandle<InverseMultimap<Index>> m_inputSimHitMeasurementsMap{
      this, "InputSimHitMeasurementsMap"};

  /// protect multi-threaded writes
  std::mutex m_writeMutex;
  /// the output file
  TFile* m_outputFile = nullptr;

  std::optional<Acts::Experimental::Histogram1> m_measurementContributingHits;
  std::optional<Acts::Experimental::Histogram1>
      m_measurementContributingParticles;
  std::optional<Acts::Experimental::Histogram1> m_measurementPurity;
  std::optional<Acts::Experimental::Histogram1> m_measurementClassification;

  std::optional<Acts::Experimental::Efficiency1> m_hitEffVsZ;
  std::optional<Acts::Experimental::Efficiency1> m_hitEffVsR;
  std::optional<Acts::Experimental::Efficiency1> m_hitEffVsEta;
  std::optional<Acts::Experimental::Efficiency1> m_hitEffVsPhi;
};

}  // namespace ActsExamples

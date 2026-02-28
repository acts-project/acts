// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/DispatchAlgorithms/DispatchEdm.hpp"
#include "ActsExamples/EventData/Cluster.hpp"
#include "ActsExamples/EventData/Measurement.hpp"
#include "ActsExamples/EventData/ProtoTrack.hpp"
#include "ActsExamples/EventData/SimParticle.hpp"
#include "ActsExamples/Framework/DataHandle.hpp"
#include "ActsExamples/Framework/IAlgorithm.hpp"
#include "ActsExamples/Framework/ProcessCode.hpp"

#include <functional>
#include <memory>
#include <string>
#include <vector>

namespace Acts {
class TrackingGeometry;
}

namespace ActsExamples {

/// This Algorithm retrieves Acts event data for pattern recognition
/// and dispatches it to a user defined function for further processing
class PatternDispatchAlgorithm final : public IAlgorithm {
 public:
  /// Configuration class it allows to connect to python functions
  class Config {
   public:
    /// @brief The pattern function that is called with the dispatch EDM and returns
    /// a vector of tracks if inputParticles and inputParticleMeasurementsMap
    /// are provided, the truth information is also passed to the pattern
    /// function, otherwise an empty structure is passed
    std::function<std::vector<DispatchTrack>(
        const DispatchMeasurements&, const DispatchParticles&,
        const DispatchParticleMeasurementsMap&)>
        patternFunction;

    /// The tracking geometry
    std::shared_ptr<const Acts::TrackingGeometry> trackingGeometry = nullptr;
    /// The input truth particles - if left free no truth association is done
    std::string inputParticles = "";
    /// The input particle-measurements map collection.
    std::string inputParticleMeasurementsMap = "";
    /// The input measurement collection
    std::string inputMeasurements = "";
    /// The output proto tracks collection.
    std::string outputProtoTracks = "";
  };

  /// Construct the smearing algorithm.
  ///
  /// @param config is the algorithm configuration
  /// @param level is the logging level
  PatternDispatchAlgorithm(Config config, Acts::Logging::Level level);

  /// Build measurement from simulation hits at input.
  ///
  /// @param ctx is the algorithm context with event information
  /// @return a process code indication success or failure
  ProcessCode execute(const AlgorithmContext& ctx) const override;

  /// Get const access to the config
  const Config& config() const { return m_cfg; }

 private:
  /// Configuration of the Algorithm
  Config m_cfg;

  ReadDataHandle<MeasurementContainer> m_inputMeasurements{this,
                                                           "InputMeasurements"};

  ReadDataHandle<SimParticleContainer> m_inputParticles{this, "InputParticles"};

  ReadDataHandle<InverseMultimap<SimBarcode>> m_inputParticleMeasurementsMap{
      this, "InputParticleMeasurementsMap"};

  WriteDataHandle<ProtoTrackContainer> m_outputProtoTracks{this,
                                                           "OutputProtoTracks"};
};

}  // namespace ActsExamples

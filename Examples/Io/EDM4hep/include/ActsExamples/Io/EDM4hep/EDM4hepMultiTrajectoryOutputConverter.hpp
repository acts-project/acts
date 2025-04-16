// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "ActsExamples/EventData/Trajectories.hpp"
#include "ActsExamples/Framework/DataHandle.hpp"
#include "ActsExamples/Io/EDM4hep/EDM4hepOutputConverter.hpp"
#include "ActsExamples/Io/Podio/CollectionBaseWriteHandle.hpp"
#include "ActsFatras/EventData/Barcode.hpp"

#include <string>

namespace ActsExamples {

/// Write out the tracks reconstructed using Combinatorial Kalman Filter to
/// EDM4hep objects.
///
/// Inpersistent information:
/// - trajectory state incomplete
/// - relation to the particles
///
/// Known issues:
/// - curvature parameter
/// - track state local coordinates are written to (D0,Z0)
/// - covariance incorrect
class EDM4hepMultiTrajectoryOutputConverter : public EDM4hepOutputConverter {
 public:
  struct Config {
    /// Input trajectory collection
    std::string inputTrajectories;
    /// Input hit-particles map collection
    std::string inputMeasurementParticlesMap;
    /// Where to place output file
    std::string outputTracks;
    /// B field in the longitudinal direction
    double Bz{};
    /// Particle hypothesis
    Acts::ParticleHypothesis particleHypothesis =
        Acts::ParticleHypothesis::pion();
  };

  /// constructor
  /// @param config is the configuration object
  /// @param level is the output logging level
  explicit EDM4hepMultiTrajectoryOutputConverter(
      const Config& config, Acts::Logging::Level level = Acts::Logging::INFO);

  /// Readonly access to the config
  const Config& config() const { return m_cfg; }

  /// Readonly access to the collections
  std::vector<std::string> collections() const final;

 protected:
  /// @brief Write method called by the base class
  /// @param [in] context is the algorithm context for consistency
  ProcessCode execute(const AlgorithmContext& context) const final;

 private:
  Config m_cfg;

  ReadDataHandle<IndexMultimap<ActsFatras::Barcode>>
      m_inputMeasurementParticlesMap{this, "InputMeasurementParticlesMaps"};

  ReadDataHandle<TrajectoriesContainer> m_inputTrajectories{
      this, "InputTrajectories"};

  CollectionBaseWriteHandle m_outputTracks{this, "OutputTracks"};
};

}  // namespace ActsExamples

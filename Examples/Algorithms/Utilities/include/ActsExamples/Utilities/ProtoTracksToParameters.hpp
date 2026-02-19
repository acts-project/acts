// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/MagneticField/MagneticFieldProvider.hpp"
#include "ActsExamples/EventData/ProtoTrack.hpp"
#include "ActsExamples/EventData/SimSeed.hpp"
#include "ActsExamples/EventData/Track.hpp"
#include "ActsExamples/Framework/DataHandle.hpp"
#include "ActsExamples/Framework/IAlgorithm.hpp"

namespace ActsExamples {

/// This algorithm computes track parameters for proto tracks.
/// Unfortunately, it is not possible to create first seeds from proto tracks,
/// and then use the more generic TrackParamEstimationAlgorithm.
///
/// The reason is, that in the end we need a consistent collection of parameters
/// and proto tracks. But if the parameter estimation fails for individual
/// seeds, it is not possible to recover a matching set of proto tracks
/// afterwards. Therefore, these two steps must be unified into one algorithm.
///
class ProtoTracksToParameters final : public IAlgorithm {
 public:
  struct Config {
    /// The proto track for that parameters should be computed
    std::string inputProtoTracks;
    /// The space point collection
    std::string inputSpacePoints;
    /// The seeds created on-the-fly from which the parameters actually are
    /// computed
    std::string outputSeeds = "seeds-from-protoTracks";
    /// The proto tracks for which parameters where computed successfully
    std::string outputProtoTracks = "remaining-protoTracks";
    /// The output parameters
    std::string outputParameters = "parameters";

    /// Whether to make tight seeds (closest 3 hits to beampipe) or large
    /// seeds (evenly spread across the proto track)
    bool buildTightSeeds = false;

    /// The tracking geometry
    std::shared_ptr<const Acts::TrackingGeometry> geometry;
    /// Magnetic field variant.
    std::shared_ptr<const Acts::MagneticFieldProvider> magneticField;

    /// The minimum magnetic field to trigger the track parameters estimation
    double bFieldMin = 0.1 * Acts::UnitConstants::T;
    /// Initial covariance matrix diagonal.
    std::array<double, 6> initialSigmas = {
        25 * Acts::UnitConstants::um,       100 * Acts::UnitConstants::um,
        0.02 * Acts::UnitConstants::degree, 0.02 * Acts::UnitConstants::degree,
        0.1 / Acts::UnitConstants::GeV,     10 * Acts::UnitConstants::ns};
    /// Inflate initial covariance.
    std::array<double, 6> initialVarInflation = {1., 1., 1., 1., 1., 1.};
    /// Particle hypothesis.
    Acts::ParticleHypothesis particleHypothesis =
        Acts::ParticleHypothesis::pion();
  };

  /// Construct the algorithm.
  ///
  /// @param cfg is the algorithm configuration
  /// @param lvl is the logging level
  ProtoTracksToParameters(Config cfg, Acts::Logging::Level lvl);

  ~ProtoTracksToParameters() override;

  /// Run the algorithm.
  ///
  /// @param ctx is the algorithm context with event information
  /// @return a process code indication success or failure
  ProcessCode execute(const AlgorithmContext& ctx) const final;

  /// Const access to the config
  const Config& config() const { return m_cfg; }

 private:
  Config m_cfg;
  Acts::BoundMatrix m_covariance = Acts::BoundMatrix::Zero();

  WriteDataHandle<SimSeedContainer> m_outputSeeds{this, "OutputSeeds"};
  WriteDataHandle<ProtoTrackContainer> m_outputProtoTracks{this,
                                                           "OutputProtoTracks"};
  WriteDataHandle<TrackParametersContainer> m_outputParameters{
      this, "OutputParameters"};
  ReadDataHandle<SimSpacePointContainer> m_inputSpacePoints{this,
                                                            "InputSpacePoints"};
  ReadDataHandle<ProtoTrackContainer> m_inputProtoTracks{this,
                                                         "InputProtoTracks"};
};

}  // namespace ActsExamples

// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/Units.hpp"
#include "Acts/EventData/ParticleHypothesis.hpp"
#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/MagneticField/MagneticFieldProvider.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/EventData/ProtoTrack.hpp"
#include "ActsExamples/EventData/Seed.hpp"
#include "ActsExamples/EventData/SeedSpacePointSelection.hpp"
#include "ActsExamples/EventData/Track.hpp"
#include "ActsExamples/Framework/DataHandle.hpp"
#include "ActsExamples/Framework/IAlgorithm.hpp"
#include "ActsExamples/Framework/ProcessCode.hpp"

#include <array>
#include <cstddef>
#include <functional>
#include <memory>
#include <string>

namespace ActsExamples {

/// Estimate track parameters from the space points of a seed.
///
/// Seeds without an estimate are dropped. Optional input seeds and proto
/// tracks are written back out filtered the same way.
class TrackParamsEstimationAlgorithm final : public IAlgorithm {
 public:
  /// Relative weight of a space point at a global position in the helix fit.
  /// Only ratios matter.
  using SpacePointWeight = std::function<double(const Acts::Vector3&)>;

  /// Weight a space point by `1 / r^exponent`.
  ///
  /// @param exponent the power of the radius to divide by
  /// @return the weight function
  static SpacePointWeight inverseRadiusPowerWeight(double exponent);

  struct Config {
    /// Input seeds collection.
    std::string inputSeeds;
    /// Input proto tracks (optional).
    std::optional<std::string> inputProtoTracks;
    /// Input particle hypothesis (optional). If not given, the static particle
    /// hypothesis from the config is used.
    std::optional<std::string> inputParticleHypotheses;
    /// Output estimated track parameters collection.
    std::string outputTrackParameters;
    /// Output seed collection - only seeds with successful parameter estimation
    /// are propagated (optional)
    std::optional<std::string> outputSeeds;
    /// Output proto track collection - only tracks with successful parameter
    /// estimation are propagated (optional)
    std::optional<std::string> outputProtoTracks;

    /// Tracking geometry for surface lookup.
    std::shared_ptr<const Acts::TrackingGeometry> trackingGeometry;
    /// Magnetic field variant.
    std::shared_ptr<const Acts::MagneticFieldProvider> magneticField;

    /// The minimum magnetic field to trigger the track parameters estimation
    double bFieldMin = 0.1 * Acts::UnitConstants::T;

    /// Which space points of the seed feed the estimate. Seeds for which the
    /// selection cannot be made are skipped.
    SeedSpacePointSelection spacePointSelection =
        SeedSpacePointSelection::FirstThree;
    /// Minimum transverse distance between the selected space points. Only the
    /// triplet selections apply it.
    double minTransverseDistance = 10 * Acts::UnitConstants::mm;
    /// Geometric refinement iterations of the circle fit. Only
    /// @c SeedSpacePointSelection::All fits a circle, a triplet is exact.
    std::size_t geometricRefineIterations = 0;
    /// Optional space point weight. Unset weights every point the same and,
    /// like the refinement, only @c SeedSpacePointSelection::All reads it.
    SpacePointWeight spacePointWeight;

    /// Initial sigmas for the track parameters.
    std::array<double, 6> initialSigmas = {
        1 * Acts::UnitConstants::mm,
        1 * Acts::UnitConstants::mm,
        1 * Acts::UnitConstants::degree,
        1 * Acts::UnitConstants::degree,
        0 * Acts::UnitConstants::e / Acts::UnitConstants::GeV,
        1 * Acts::UnitConstants::ns};
    /// Initial sigma(q/pt) for the track parameters.
    /// @note The resulting q/p sigma is added to the one in `initialSigmas`
    double initialSigmaQoverPt =
        0.1 * Acts::UnitConstants::e / Acts::UnitConstants::GeV;
    /// Initial sigma(pt)/pt for the track parameters.
    /// @note The resulting q/p sigma is added to the one in `initialSigmas`
    double initialSigmaPtRel = 0.1;
    /// Inflate initial covariance.
    std::array<double, 6> initialVarInflation = {1., 1., 1., 1., 1., 1.};
    /// Inflate time covariance if no time measurement is available.
    double noTimeVarInflation = 100.;

    /// Particle hypothesis.
    Acts::ParticleHypothesis particleHypothesis =
        Acts::ParticleHypothesis::pion();
  };

  /// Construct the track parameters making algorithm.
  ///
  /// @param cfg is the algorithm configuration
  /// @param lvl is the logging level
  explicit TrackParamsEstimationAlgorithm(
      const Config& cfg, std::unique_ptr<const Acts::Logger> logger = nullptr);

  /// Run the track parameters making algorithm.
  ///
  /// @param ctx is the algorithm context with event information
  /// @return a process code indication success or failure
  ProcessCode execute(const AlgorithmContext& ctx) const override;

  /// Const access to the config
  const Config& config() const { return m_cfg; }

 private:
  Config m_cfg;

  ReadDataHandle<SeedContainer> m_inputSeeds{this, "InputSeeds"};
  ReadDataHandle<ProtoTrackContainer> m_inputTracks{this, "InputTracks"};
  ReadDataHandle<std::vector<Acts::ParticleHypothesis>>
      m_inputParticleHypotheses{this, "InputParticleHypotheses"};

  WriteDataHandle<TrackParametersContainer> m_outputTrackParameters{
      this, "OutputTrackParameters"};
  WriteDataHandle<SeedContainer> m_outputSeeds{this, "OutputSeeds"};
  WriteDataHandle<ProtoTrackContainer> m_outputTracks{this, "OutputTracks"};
};

}  // namespace ActsExamples

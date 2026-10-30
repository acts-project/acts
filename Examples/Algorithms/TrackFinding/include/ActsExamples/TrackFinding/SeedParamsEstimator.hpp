// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/Definitions/Units.hpp"
#include "Acts/EventData/ParticleHypothesis.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/MagneticField/MagneticFieldContext.hpp"
#include "Acts/MagneticField/MagneticFieldProvider.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "Acts/Utilities/Result.hpp"
#include "ActsExamples/EventData/Seed.hpp"
#include "ActsExamples/EventData/SeedSpacePointSelection.hpp"
#include "ActsExamples/EventData/SpacePoint.hpp"
#include "ActsExamples/EventData/Track.hpp"

#include <array>
#include <cstddef>
#include <functional>
#include <memory>
#include <vector>

namespace ActsExamples {

/// Estimate track parameters from the space points of a seed, on the surface
/// of its first or its last space point. The fit and the space point selection
/// start at that space point.
class SeedParamsEstimator {
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

  /// The space point of the seed that the estimate is given on
  enum class Anchor { First, Last };

  /// Why a seed has no estimate
  enum class Failure {
    SpacePointSelection,
    MissingSurface,
    FieldLookup,
    WeakField,
    Fit,
    Degenerate,
    Transport,
  };

  /// Buffers reused across seeds
  struct Cache {
    /// The magnetic field cache
    Acts::MagneticFieldProvider::Cache field;
    /// The space points of the seed from the outer end
    std::vector<SpacePointIndex> candidates{};
    /// The selected triplet
    std::array<SpacePointIndex, 3> triplet{};
    /// The positions of the selected space points
    std::vector<Acts::Vector3> positions{};
    /// The weights of the selected space points
    std::vector<double> weights{};
  };

  /// @param cfg the configuration
  /// @param logger the logger
  explicit SeedParamsEstimator(
      Config cfg, std::unique_ptr<const Acts::Logger> logger = nullptr);

  /// @param mctx the magnetic field context
  /// @return a cache for @ref estimate
  Cache makeCache(const Acts::MagneticFieldContext& mctx) const;

  /// Estimate the parameters of one seed.
  ///
  /// @param gctx the geometry context
  /// @param mctx the magnetic field context
  /// @param cache the buffers
  /// @param seed the seed
  /// @param hypothesis the particle hypothesis
  /// @param freeParams given free parameters which replace the fit, unless
  ///        they are null or have no momentum
  /// @param anchor the space point the estimate is given on
  /// @return the parameters, or why there are none
  Acts::Result<TrackParameters, Failure> estimate(
      const Acts::GeometryContext& gctx, const Acts::MagneticFieldContext& mctx,
      Cache& cache, const ConstSeedProxy& seed,
      const Acts::ParticleHypothesis& hypothesis,
      const Acts::FreeVector* freeParams, Anchor anchor) const;

  /// Get the configuration
  /// @return the configuration
  const Config& config() const { return m_cfg; }

 private:
  Config m_cfg;
  std::unique_ptr<const Acts::Logger> m_logger;
  std::shared_ptr<const Acts::Logger> m_propagatorLogger;

  const Acts::Logger& logger() const { return *m_logger; }
};

}  // namespace ActsExamples

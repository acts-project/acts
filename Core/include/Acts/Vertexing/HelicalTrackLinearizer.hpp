// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/MagneticField/MagneticFieldContext.hpp"
#include "Acts/MagneticField/MagneticFieldProvider.hpp"
#include "Acts/MagneticField/NullBField.hpp"
#include "Acts/Propagator/EigenStepper.hpp"
#include "Acts/Propagator/Propagator.hpp"
#include "Acts/Utilities/Result.hpp"
#include "Acts/Vertexing/LinearizedTrack.hpp"

namespace Acts {

/// @class HelicalTrackLinearizer
/// Linearizes the track parameters at the PCA to a user-provided
/// point (linPoint). The track parameters are written as a function
/// of the global PCA position and the momentum of the particle at
/// the PCA. The linearization then reads (see Eq. 5.7 in Ref(1)):
///
/// q = A (r - r_0) + B (p - p_0) + c,
///
/// where q are the Perigee parameters wrt linPoint, {r_0} r is the {initial}
/// 4D PCA position, {p_0} p is the {initial} momentum (phi, theta, q/p) at the
/// PCA, and c is the constant term of the expansion. A and B are matrices of
/// derivatives, denoted hereafter as "positionJacobian" and
/// "momentumJacobian" respectively.
///
/// This class computes A and B using the analytic formulae of Ref. (1).
///
/// Ref.(1) - CERN-THESIS-2010-027, Giacinto Piacquadio (Freiburg U.)
///
/// @tparam propagator_t Propagator type
/// @tparam propagator_options_t Propagator options type
template <typename propagator_t,
          typename propagator_options_t = PropagatorOptions<>>
class HelicalTrackLinearizer {
 public:
  using Propagator_t = propagator_t;

  /// State struct
  struct State {
    /// @brief The state constructor
    ///
    /// @param fieldCacheIn The magnetic field cache
    State(MagneticFieldProvider::Cache fieldCacheIn)
        : fieldCache(std::move(fieldCacheIn)) {}
    /// Magnetic field cache
    MagneticFieldProvider::Cache fieldCache;
  };

  /// @brief Configuration struct
  struct Config {
    /// @ Config constructor if magnetic field is present
    ///
    /// @param bIn The magnetic field
    /// @param prop The propagator
    Config(std::shared_ptr<const MagneticFieldProvider> bIn,
           std::shared_ptr<const Propagator_t> prop)
        : bField(std::move(bIn)), propagator(std::move(prop)) {}

    /// @brief Config constructor without B field -> uses NullBField
    ///
    /// @param prop The propagator
    Config(std::shared_ptr<const Propagator_t> prop)
        : bField{std::make_shared<NullBField>()}, propagator(std::move(prop)) {}

    // The magnetic field
    std::shared_ptr<const MagneticFieldProvider> bField;
    // The propagator
    std::shared_ptr<const Propagator_t> propagator;

    // Minimum q/p value
    double minQoP = 1e-15;
    // Maximum curvature value
    double maxRho = 1e+15;
  };

  /// @brief Constructor
  ///
  /// @param config Configuration object
  /// @param _logger a logger instance
  HelicalTrackLinearizer(const Config& config,
                         std::unique_ptr<const Logger> _logger =
                             getDefaultLogger("HelTrkLinProp", Logging::INFO))
      : m_cfg(config), m_logger{std::move(_logger)} {}

  /// @brief Function that linearizes BoundTrackParameters at
  /// the PCA to a given Perigee surface
  ///
  /// @param params Parameters to linearize
  /// @param linPoint Point which defines the Perigee.
  /// @note Transverse plane of the Perigee corresponding to @p linPoint is
  /// parallel to the global x-y plane
  /// @param gctx Geometry context
  /// @param mctx Magnetic field context
  /// @param state Linearizer state object
  ///
  /// @return Linearized track
  Result<LinearizedTrack> linearizeTrack(const BoundTrackParameters& params,
                                         const Vector4& linPoint,
                                         const Acts::GeometryContext& gctx,
                                         const Acts::MagneticFieldContext& mctx,
                                         State& state) const;

 private:
  /// Configuration object
  const Config m_cfg;

  std::unique_ptr<const Logger> m_logger;

  const Logger& logger() const { return *m_logger; }
};

}  // namespace Acts

#include "HelicalTrackLinearizer.ipp"

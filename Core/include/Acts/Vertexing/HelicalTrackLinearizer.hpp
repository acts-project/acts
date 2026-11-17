// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Propagator/EigenStepper.hpp"
#include "Acts/Propagator/Propagator.hpp"
#include "Acts/Utilities/Delegate.hpp"
#include "Acts/Utilities/Result.hpp"
#include "Acts/Vertexing/LinearizedTrack.hpp"

#include <memory>

namespace Acts {

class MagneticFieldQuery;

/// @class HelicalTrackLinearizer
/// Linearizes the track parameters at the PCA to a user-provided
/// point (linPoint). The track parameters are written as a function
/// of the global PCA position and the momentum of the particle at
/// the PCA. The linearization then reads (see Eq. 5.7 in Ref. (1)):
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
/// Ref. (1) - CERN-THESIS-2010-027, Giacinto Piacquadio (Freiburg U.)
class HelicalTrackLinearizer {
 public:
  /// @brief Configuration struct
  struct Config {
    /// Track propagator for linearization calculations
    std::shared_ptr<const BasePropagator> propagator;

    /// Tolerance determining how close we need to get to the Perigee surface to
    /// reach it during propagation
    double targetTolerance = 1e-12;
  };

  /// @brief Constructor
  ///
  /// @param config Configuration object
  /// @param _logger a logger instance
  explicit HelicalTrackLinearizer(
      const Config& config,
      std::unique_ptr<const Logger> _logger = getDefaultLogger("HelTrkLinProp",
                                                               Logging::INFO))
      : m_cfg(config), m_logger{std::move(_logger)} {
    if (!m_cfg.propagator) {
      throw std::invalid_argument("HelicalTrackLinearizer: propagator is null");
    }
  }

  /// Linearize using a caller-owned field query for the PCA lookup.
  /// The propagator keeps its own cache and must use the same physical field.
  /// @param params Parameters to linearize
  /// @param linPointTime Time of the linearization point
  /// @param perigeeSurface Surface defining the linearization point
  /// @note The transverse plane of @p perigeeSurface is parallel to the
  /// global x-y plane
  /// @param gctx Geometry context
  /// @param fieldQuery Query providing both the PCA field and propagation context
  /// @return Linearized track or a propagation/field error
  Result<LinearizedTrack> linearizeTrack(const BoundTrackParameters& params,
                                         double linPointTime,
                                         const Surface& perigeeSurface,
                                         const GeometryContext& gctx,
                                         MagneticFieldQuery& fieldQuery) const;

 private:
  /// Configuration object
  const Config m_cfg;

  std::unique_ptr<const Logger> m_logger;

  const Logger& logger() const { return *m_logger; }
};

}  // namespace Acts

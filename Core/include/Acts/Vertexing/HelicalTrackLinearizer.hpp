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
/// Linearizes the measurement equation (dependance of track
/// parameters on the vertex position and track momentum at vertex)
/// at the vicinity of the user-provided linearization point.
///
/// The measurement equation is linearized in the following way:
///
/// q_k= A_k (x_k - x_0k) + B_k (p_k - p_0k) + c_k
///
/// where q_k are the parameters at perigee nearest to the lin point,
/// x_k is the position of the vertex, p_k the track momentum at the vertex,
/// and c_k is the constant term of expansion. A_k and B_k are matrices
/// of derivatives, denoted hereafter as "positionJacobian" and
/// "momentumJacobian" respectively.
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
  /// given linearization point
  ///
  /// @param params Parameters to linearize
  /// @param linPoint Linearization point
  /// @param gctx The geometry context
  /// @param mctx The magnetic field context
  /// @param state The state object
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

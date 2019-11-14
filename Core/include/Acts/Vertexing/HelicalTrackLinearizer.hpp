// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/MagneticField/MagneticFieldContext.hpp"
#include "Acts/Propagator/EigenStepper.hpp"
#include "Acts/Propagator/Propagator.hpp"
#include "Acts/Utilities/Definitions.hpp"
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
  using PropagatorOptions_t = propagator_options_t;

 public:
  using Propagator_t = propagator_t;
  using BField_t = typename Propagator_t::Stepper::BField;

  /// @brief Helper function to set up the correct PropagatorOptions with
  /// propagation direction set to backward. To be used when setting up the
  /// propagator options for the HelicalTrackLinearizer Config.
  ///
  /// @param gc The GeometryContext
  /// @param mv The MagneticFieldContext
  ///
  /// @return The PropagatorOptions with direction = backward
  static PropagatorOptions_t getDefaultPropagatorOptions(
      const GeometryContext& gc, const MagneticFieldContext& mc) {
    PropagatorOptions_t options(gc, mc);
    options.direction = backward;
    return options;
  }

  /// @brief Configuration struct
  ///
  /// @param bIn The magnetic field
  /// @param prop The propagator
  /// @param propOptions The propagator options
  struct Config {
    Config(const BField_t& bIn, std::shared_ptr<Propagator_t> prop,
           PropagatorOptions_t propOptions)
        : bField(bIn), propagator(std::move(prop)), pOptions(propOptions) {
      assert(pOptions.direction == backward);
    }

    /// @brief Config constructor if BField_t == int (no B-Field provided),
    ///        sets int bField to 0
    template <typename T = BField_t,
              std::enable_if_t<std::is_same<T, int>::value, int> = 0>
    Config(std::shared_ptr<Propagator_t> prop, PropagatorOptions_t propOptions)
        : bField(0), propagator(std::move(prop)), pOptions(propOptions) {
      assert(pOptions.direction == backward);
    }

    // The magnetic field
    BField_t bField;
    // The propagator
    std::shared_ptr<Propagator_t> propagator;
    // The propagator options
    PropagatorOptions_t pOptions;
  };

  /// @brief Constructor
  ///
  /// @param config Configuration object
  HelicalTrackLinearizer(const Config& config) : m_cfg(config) {}

  /// @brief Function that linearizes BoundParameters at
  /// given linearization point
  ///
  /// @param params Parameters to linearize
  /// @param linPoint Linearization point
  ///
  /// @return Linearized track
  Result<LinearizedTrack> linearizeTrack(
      const BoundParameters* params, const SpacePointVector& linPoint) const;

 private:
  /// @brief Method that returns the magnetic field value at a given position
  ///        Enabled if BField_t == int (no B-Field provided), returns 0.
  ///
  /// @param linPointPos Position for which to get the magnetic field value
  /// @return The magnetic field value
  template <typename T = BField_t,
            std::enable_if_t<std::is_same<T, int>::value, int> = 0>
  double getBField(const Acts::Vector3D& /*linPointPos*/) const {
    return m_cfg.bField;
  }

  /// @brief Method that returns the magnetic field value at a given position
  ///        Enabled if BField_t != int (B-Field provided)
  ///
  /// @param linPointPos Position for which to get the magnetic field value
  /// @return The magnetic field value
  template <typename T = BField_t,
            std::enable_if_t<!std::is_same<T, int>::value, int> = 0>
  double getBField(const Acts::Vector3D& linPointPos) const {
    return m_cfg.bField.getField(linPointPos)[eZ];
  }

  /// Configuration object
  const Config m_cfg;
};

}  // namespace Acts

#include "HelicalTrackLinearizer.ipp"

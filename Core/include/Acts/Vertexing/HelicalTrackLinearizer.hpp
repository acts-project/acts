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
/// @tparam bfield_t Magnetic field type
/// @tparam propagator_t Propagator type
/// @tparam action_list_t Propagator action list type
/// @tparam aborter_list_t Propagator aborter list type
template <typename bfield_t,
          typename propagator_t = Propagator<EigenStepper<bfield_t>>,
          typename action_list_t = ActionList<>,
          typename aborter_list_t = AbortList<>>
class HelicalTrackLinearizer {
 public:
  using Propagator_t = propagator_t;

  /// @brief Helper function to set up the correct PropagatorOptions with
  /// propagation direction set to backward. To be used when setting up the
  /// propagator options for the HelicalTrackLinearizer Config.
  ///
  /// @param gc The GeometryContext
  /// @param mv The MagneticFieldContext
  ///
  /// @return The PropagatorOptions with direction = backward
  static PropagatorOptions<action_list_t, aborter_list_t>
  getDefaultPropagatorOptions(const GeometryContext& gc,
                              const MagneticFieldContext& mc) {
    PropagatorOptions<action_list_t, aborter_list_t> options(gc, mc);
    options.direction = backward;
    return options;
  }

  struct Config {
    Config(const bfield_t& bIn, const Propagator_t& prop,
           PropagatorOptions<action_list_t, aborter_list_t> propOptions)
        : bField(bIn), propagator(prop), pOptions(propOptions) {
      assert(pOptions.direction == backward);
    }

    bfield_t bField;

    Propagator_t propagator;

    PropagatorOptions<action_list_t, aborter_list_t> pOptions;
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
  /// Configuration object
  const Config m_cfg;
};

}  // namespace Acts

#include "HelicalTrackLinearizer.ipp"

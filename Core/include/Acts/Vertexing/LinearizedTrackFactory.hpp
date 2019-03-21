// This file is part of the Acts project.
//
// Copyright (C) 2019 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/Propagator/Propagator.hpp"
#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Vertexing/LinearizedTrack.hpp"

namespace Acts {

/// @class LinearizedTrackFactory
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
/// @tparam BField Magnetic field type
/// @tparam Propagator_t Propagator type
/// @tparam action_list_t Propagator action list type
/// @tparam aborter_list_t Propagator aborter list type

template <typename BField,
          typename Propagator_t,
          typename action_list_t  = ActionList<>,
          typename aborter_list_t = AbortList<>>
class LinearizedTrackFactory
{
public:
  struct Config
  {
    BField bField;

    /// Default propagator options
    PropagatorOptions<action_list_t, aborter_list_t> propagatorOptions;

    Config(BField bIn) : bField(std::move(bIn)){};
  };

  /// @brief Constructor with BField
  ///
  /// @param config Configuration object
  LinearizedTrackFactory(const Config& config) : m_cfg(config) {}

  /// @brief Function that linearizes BoundParameters at
  /// given linearization point
  ///
  /// @param params Parameters to linearize
  /// @param linPoint Linearization point
  /// @param propagator Propagator
  ///
  /// @return Linearized track
  LinearizedTrack
  linearizeTrack(const BoundParameters* params,
                 const Vector3D&        linPoint,
                 const Propagator_t&    propagator) const;

private:
  // Configuration object
  Config m_cfg;
};

}  // namespace Acts

#include "LinearizedTrackFactory.ipp"

// This file is part of the Acts project.
//
// Copyright (C) 2016-2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Vertexing/LinearizedTrack.hpp"

namespace Acts {

/**
 * @class LinearizedTrackFactory
 * From ATHENA:
 * Linearizes the measurement equation (dependance of track
 * parameters on the vertex position and track mometum at vertex)
 * at the vicinity of the user-provided linearization point.
 *
 * The measurement equation is linearized in the following way:
 *
 * q_k= A_k (x_k - x_0k) + B_k (p_k - p_0k) + c_k
 *
 * where q_k are the parameters at perigee nearest to the lin point,
 * x_k is the position of the vertex, p_k the track momentum at the vertex,
 * and c_k is the constant term of expansion. A_k and B_k are matrices
 * of derivatives, denoted hereafter as "positionJacobian" and
 * "momentumJacobian"
 * respectively.
 *
 */

template <typename BField>
class LinearizedTrackFactory
{
public:
  struct Config
  {
    BField bField;

    Config(BField bIn) : bField(std::move(bIn)){};
  };

  /// Constructor with BField
  LinearizedTrackFactory(const Config& config) : m_cfg(config) {}

  /// Default destructor
  ~LinearizedTrackFactory() = default;

  LinearizedTrack*
  linearizeTrack(const Acts::BoundParameters* params,
                 const Acts::Vector3D&        linPoint) const;

private:
  // Configuration object
  Config m_cfg;
};

}  // namespace Acts

#include "LinearizedTrackFactory.ipp"

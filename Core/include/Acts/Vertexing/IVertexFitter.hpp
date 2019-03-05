// This file is part of the Acts project.
//
// Copyright (C) 2019 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Vertexing/Vertex.hpp"

namespace Acts {

/// @class IVertexFitter
///
/// @brief Virtual base class for VertexFitters
///
/// @tparam InputTrack Track object type
/// @tparam Propagator_t Propagator type

template <typename InputTrack, typename Propagator_t>
class IVertexFitter
{
public:
  /// @brief Default virtual destructor
  virtual ~IVertexFitter() = default;

  /// @param paramVector Vector of track objects to fit vertex to
  /// @param propagator Propagator
  /// @param constraint Constraint of the fit, position of constraint is
  /// starting point
  ///
  /// @return Fitted vertex
  virtual Vertex<InputTrack>
  fit(const std::vector<InputTrack>& paramVector,
      const Propagator_t&            propagator,
      Vertex<InputTrack>             constraint
      = Vertex<InputTrack>(Vector3D(0., 0., 0.))) const = 0;
};

}  // namespace Acts

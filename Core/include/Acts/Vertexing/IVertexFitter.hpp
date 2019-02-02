// This file is part of the Acts project.
//
// Copyright (C) 2016-2019 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Vertexing/Vertex.hpp"

namespace Acts {

/// @class IVertexFitter
///
/// @brief Virtual base class for VertexFitters

template <typename InputTrack>
class IVertexFitter
{
public:
  /// @param paramVector Vector of track objects to fit vertex to
  /// @param startingPoint Constraint of the fit, position of constraint is
  /// starting point
  ///
  /// @return Fitted vertex
  virtual Vertex<InputTrack>
  fit(const std::vector<InputTrack>& paramVector,
      Vertex<InputTrack>             constraint) const = 0;
};

}  // namespace Acts
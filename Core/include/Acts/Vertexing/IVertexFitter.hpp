// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/MagneticField/MagneticFieldContext.hpp"
#include "Acts/Utilities/Result.hpp"
#include "Acts/Vertexing/Vertex.hpp"
#include "Acts/Vertexing/VertexFitInput.hpp"

namespace Acts {

/// Common interface for fitting one vertex.
/// Implementations create their working state within each call. Coupled
/// multi-vertex and cascade fitting are outside this interface's contract.
class IVertexFitter {
 public:
  virtual ~IVertexFitter() = default;

  /// Fit one vertex from tracks, a seed and an optional vertex constraint.
  /// @param input Fit inputs, including the initial linearization point
  /// @param gctx Geometry context for this fit
  /// @param mctx Magnetic field context for this fit
  /// @return Fitted vertex with its fitted tracks, or a fitting error
  virtual Result<Vertex> fit(const VertexFitInput& input,
                             const GeometryContext& gctx,
                             const MagneticFieldContext& mctx) const = 0;
};

}  // namespace Acts

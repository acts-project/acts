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

class MagneticFieldQuery;

/// Common interface for fitting one vertex.
/// Callers can reuse a field query across fits. Coupled
/// multi-vertex and cascade fitting are outside this interface's contract.
class IVertexFitter {
 public:
  virtual ~IVertexFitter() = default;

  /// Fit one vertex from tracks, a seed and an optional vertex constraint.
  /// @param input Fit inputs, including the initial linearization point
  /// @param gctx Geometry context for this fit
  /// @param mctx Magnetic field context for this fit
  /// @param field Reusable field query, created for mctx; not retained by the fit
  /// @return Fitted vertex with its fitted tracks, or a fitting error
  virtual Result<Vertex> fit(const VertexFitInput& input,
                             const GeometryContext& gctx,
                             const MagneticFieldContext& mctx,
                             MagneticFieldQuery& field) const = 0;
};

}  // namespace Acts

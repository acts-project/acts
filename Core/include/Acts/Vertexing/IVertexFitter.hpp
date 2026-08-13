// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/MagneticField/MagneticFieldContext.hpp"
#include "Acts/Utilities/Any.hpp"
#include "Acts/Utilities/Result.hpp"
#include "Acts/Vertexing/VertexFitProblem.hpp"
#include "Acts/Vertexing/VertexingOptions.hpp"

#include <vector>

namespace Acts {

class Vertex;
struct InputTrack;
struct VertexingOptions;

/// @brief Common interface for vertex fitters
///
/// A fitter takes a @c VertexFitProblem — the vertices to fit, their
/// constraints and the tracks assigned to them — and updates the vertices in
/// place. Scratch data is kept in a type-erased, caller-owned cache, so that a
/// fitter instance is stateless and can be shared between threads.
///
/// @note Implementations differ in one semantic property that this interface
/// deliberately does not hide: **whether the fitter couples vertices that
/// share tracks**. @c AdaptiveMultiVertexFitter does — it fits all vertices in
/// the problem simultaneously and competes shared tracks between them.
/// @c FullBilloirVertexFitter does not — it fits each vertex independently.
/// Substituting one for the other is therefore a change of algorithm, not of
/// configuration. For this reason @c AdaptiveMultiVertexFinder stays bound to
/// the concrete @c AdaptiveMultiVertexFitter rather than taking this interface.
class IVertexFitter {
 public:
  /// Move-only type-erased wrapper for concrete fitter cache objects
  using Cache = Acts::AnyBase<128, false>;

  /// Function to create a cache object for this concrete vertex fitter
  /// @param mctx The magnetic field context
  /// @return The cache object
  virtual Cache makeCache(const MagneticFieldContext& mctx) const = 0;

  /// Fit all vertices in @p problem, updating them in place
  /// @param problem The multi-vertex fit problem
  /// @param vertexingOptions The vertexing options
  /// @param cache The cache object (needs to be created via @c makeCache)
  /// @return Result indicating success or failure of the fit
  virtual Result<void> fit(VertexFitProblem& problem,
                           const VertexingOptions& vertexingOptions,
                           Cache& cache) const = 0;

  /// Add vertices to an existing fit problem and fit them
  ///
  /// The new vertices must already have their candidate entries filled in
  /// @p problem. Fitters that couple vertices sharing tracks will also refit
  /// the vertices connected to the new ones.
  ///
  /// @param problem The multi-vertex fit problem
  /// @param newVertices The vertices to add to the fit
  /// @param vertexingOptions The vertexing options
  /// @param cache The cache object (needs to be created via @c makeCache)
  /// @return Result indicating success or failure of the fit
  virtual Result<void> addVertices(VertexFitProblem& problem,
                                   const std::vector<Vertex*>& newVertices,
                                   const VertexingOptions& vertexingOptions,
                                   Cache& cache) const = 0;

  /// Convenience method for the common single-vertex case
  ///
  /// Fits one vertex to @p trackVector without a constraint and returns it.
  ///
  /// @param trackVector The tracks to fit a vertex to
  /// @param vertexingOptions The vertexing options
  /// @param cache The cache object (needs to be created via @c makeCache)
  /// @return The fitted vertex
  virtual Result<Vertex> fitSingle(const std::vector<InputTrack>& trackVector,
                                   const VertexingOptions& vertexingOptions,
                                   Cache& cache) const = 0;

  /// Virtual destructor
  virtual ~IVertexFitter() = default;
};

}  // namespace Acts

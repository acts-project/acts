// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "Acts/Material/ISurfaceMaterial.hpp"
#include "Acts/Material/MaterialInteraction.hpp"
#include "Acts/Material/TrackingGeometryMaterial.hpp"
#include "Acts/Material/interface/IAssignmentFinder.hpp"

#include <map>
#include <memory>
#include <span>
#include <vector>

namespace Acts {

class Surface;
namespace detail {
class WorkerPool;
}  // namespace detail

/// @brief Interface for the material mapping, this is the accumulation step
/// @ingroup material_mapping
class ISurfaceMaterialAccumulator {
 public:
  /// The state of the material accumulator, this is used
  /// to cache information across tracks/events
  class State {
   public:
    virtual ~State() = default;
  };

  /// Virtual destructor
  virtual ~ISurfaceMaterialAccumulator() = default;

  /// Factory for creating the state
  /// @param gctx the geometry context
  /// @return Unique pointer to a new state object for material accumulation
  virtual std::unique_ptr<State> createState(
      const GeometryContext& gctx) const = 0;

  /// @brief Accumulate the material interaction on the surface
  ///
  /// @param state is the state of the accumulator
  /// @param gctx the geometry context
  /// @param interactions is the material interactions, with assigned surfaces
  /// @param surfacesWithoutAssignment are the surfaces without assignment
  ///
  /// @note this the track average over the binned material
  virtual void accumulate(
      State& state, const GeometryContext& gctx,
      const std::vector<MaterialInteraction>& interactions,
      const std::vector<IAssignmentFinder::SurfaceAssignment>&
          surfacesWithoutAssignment) const = 0;

  /// The accumulation input of a single track
  struct TrackAccumulation {
    /// The material interactions, with assigned surfaces
    const std::vector<MaterialInteraction>* interactions = nullptr;
    /// The candidate surfaces without assignment
    const std::vector<IAssignmentFinder::SurfaceAssignment>*
        surfacesWithoutAssignment = nullptr;
  };

  /// @brief Accumulate the material of a batch of tracks
  ///
  /// The result must be identical to calling the single-track
  /// @ref accumulate for every track in order. The default does exactly that
  /// and ignores @p pool; implementations can override this to process
  /// independent parts of the accumulation concurrently.
  ///
  /// @param state is the state of the accumulator
  /// @param gctx the geometry context
  /// @param tracks the per-track accumulation inputs, in mapping order
  /// @param pool optional worker pool, nullptr means serial execution
  virtual void accumulate(State& state, const GeometryContext& gctx,
                          std::span<const TrackAccumulation> tracks,
                          detail::WorkerPool* pool) const {
    static_cast<void>(pool);
    for (const auto& track : tracks) {
      accumulate(state, gctx, *track.interactions,
                 *track.surfacesWithoutAssignment);
    }
  }

  /// Finalize the surface material maps
  ///
  /// @param state the state of the accumulator
  /// @param gctx the geometry context
  ///
  /// @note this does the run average over the (binned) material
  /// @return Map of geometry IDs to finalized surface material objects
  virtual std::map<GeometryIdentifier, std::shared_ptr<const ISurfaceMaterial>>
  finalizeMaterial(State& state, const GeometryContext& gctx) const = 0;

  /// Finalize material assignments, including stable keys when supported.
  /// The default adapts the existing ID-only finalizeMaterial implementation.
  /// @param state State created by this accumulator
  /// @param gctx Geometry context
  /// @return Finalized material maps
  virtual TrackingGeometryMaterial finalizeMaps(
      State& state, const GeometryContext& gctx) const {
    return {finalizeMaterial(state, gctx), {}};
  }
};

}  // namespace Acts

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
#include "Acts/Material/interface/IAssignmentFinder.hpp"

#include <map>
#include <memory>
#include <vector>

namespace Acts {

class Surface;

/// @brief Interface for the material mapping, this is the accumulation step
class ISurfaceMaterialAccumulater {
 public:
  /// The state of the material accumulater, this is used
  /// to cache information across tracks/events
  class State {
   public:
    virtual ~State() = default;
  };

  /// Virtual destructor
  virtual ~ISurfaceMaterialAccumulater() = default;

  /// Factory for creating the state
  /// @param gctx the geometry context
  /// @return Unique pointer to a new state object for material accumulation
  virtual std::unique_ptr<State> createState(
      const GeometryContext& gctx) const = 0;

  /// @brief Accumulate the material interaction on the surface
  ///
  /// @param state is the state of the accumulater
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

  /// Finalize the surface material maps
  ///
  /// @param state the state of the accumulator
  /// @param gctx the geometry context
  ///
  /// @note this does the run average over the (binned) material
  /// @return Map of geometry IDs to finalized surface material objects
  virtual std::map<GeometryIdentifier, std::shared_ptr<const ISurfaceMaterial>>
  finalizeMaterial(State& state, const GeometryContext& gctx) const = 0;
};

}  // namespace Acts

// This file is part of the Acts project.
//
// Copyright (C) 2024 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/GeometryHierarchyMap.hpp"
#include "Acts/Material/MaterialInteraction.hpp"

#include <array>
#include <optional>
#include <tuple>
#include <vector>

namespace Acts {

class Surface;

namespace MaterialInteractionAssignment {

/// Surface assignement definition
using SurfaceAssignment = std::tuple<const Surface*, Vector3, Vector3>;

/// @brief definition of a global veo on a material interaction
using GlobalVeto = std::function<bool(const MaterialInteraction&)>;

/// @brief definition of a local veto on a material interaction
///
/// This can take already the suggested surface assignment into account
/// return true if the assignment should be vetoed
using LocalVeto =
    std::function<bool(const MaterialInteraction&, const SurfaceAssignment&)>;

/// @brief definition of possible re-assignments to next, this could e.g.
/// be used for respecting pre/post mappint directives that are not fully
/// handled by closest distance matching
///
/// The provided parameters are the mutable material interaction, the suggested
/// assignment and the next possible assignment, due to the ordered nature of
/// the material interactions, assignment to previous is excluded
///
/// @note this changes the MaterialInteraction if the re-assignment is accepted
using ReAssignment = std::function<void(
    MaterialInteraction&, const SurfaceAssignment&, const SurfaceAssignment&)>;

/// @brief Options for the material interaction matcher
/// The options are used to specify the vetos for the assignment
struct Options {
  /// Allow global vetos for the assignment, e.g. restricting
  /// the assignment to a specific volume
  std::vector<GlobalVeto> globalVetos = {};

  /// Allow specific vetoes for the assignment, e.g. only locally to
  /// surface, rest to next mapper (e.g. volume mapper)
  GeometryHierarchyMap<LocalVeto> localVetos = {};

  /// Allow re-assignment of the material interaction, e.g. to respect
  GeometryHierarchyMap<ReAssignment> reAssignments = {};
};

/// @brief Match the material interactions to surfaces intersections while respecting
/// eventual vetos for the assignment
///
/// @param gctx is the geometry context
/// @param materialInteractions is the vector of material interaction
/// @param intersectedSurfaces is the geometry identifier
/// @param options are the options for the assignment
///
/// @return a tuple of [ Assigned, NotAssigned, Surfaces with no assignments ]
std::tuple<std::vector<MaterialInteraction>, std::vector<MaterialInteraction>,
           std::vector<SurfaceAssignment>>
assign(const GeometryContext& gctx,
       const std::vector<MaterialInteraction>& materialInteractions,
       const std::vector<SurfaceAssignment>& intersectedSurfaces,
       const Options& options = Options());

}  // namespace MaterialInteractionAssignment

}  // namespace Acts

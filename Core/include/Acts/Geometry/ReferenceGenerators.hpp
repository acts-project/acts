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
#include "Acts/Geometry/IReferenceGenerator.hpp"
#include "Acts/Geometry/Polyhedron.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/AxisDefinitions.hpp"

#include <ranges>
#include <vector>

namespace Acts {

/// A struct to access the center position as a sole reference
///
/// This generator will provide only one filling point and hence
/// only a single bin in the indexed grid.
struct CenterReferenceGenerator : public IReferenceGenerator {
  /// Helper to access the Center point of for filling the grid
  ///
  /// @param gctx the geometry context of this operation
  /// @param surface the surface for which the reference point is to be accessed
  ///
  /// @return a vector of reference points for filling
  const std::vector<Vector3> references(const GeometryContext& gctx,
                                        const Surface& surface) const override {
    return {surface.center(gctx)};
  }
};

/// A struct to access reference positions based on bin values
///
/// @tparam bVAL the binning value to be used for the binning position call
///
/// This generator will provide only one filling point and hence
/// only a single bin in the indexed grid.
template <AxisDirection bVAL>
struct AxisDirectionReferenceGenerator : public IReferenceGenerator {
  /// Helper to access a reference position based on binning value
  ///
  /// @param gctx the geometry context of this operation
  /// @param surface the surface for which the reference point is to be accessed
  ///
  /// @return a vector of reference points for filling
  const std::vector<Vector3> references(const GeometryContext& gctx,
                                        const Surface& surface) const override {
    return {surface.referencePosition(gctx, bVAL)};
  }
};

/// A struct to access generated vertices from surface polyhedrons
/// These vertices are then used to find the bin boundary box for the
/// indexed grid.
///
///
/// The grid filling then completes the empty bins in between and
/// expands if necessary.
struct PolyhedronReferenceGenerator : public IReferenceGenerator {
  /// This is for the barycenter addition
  bool addBarycenter = false;

  /// @brief  The number of segments for the polyhedron approximation
  int nSegements = 1;

  /// Absolute expansion value for the reference points
  double expansionValue = 0.0;

  /// Helper to access the Center point of for filling the grid
  ///
  /// @param gctx the geometry context of this operation
  /// @param surface the surface for which the reference point is to be accessed
  ///
  /// @return a vector of reference points for filling
  const std::vector<Vector3> references(const GeometryContext& gctx,
                                        const Surface& surface) const override;
};

/// A Projected reference generator which projects the polyhedron vertices onto
/// a given reference surface.
///
struct ProjectedReferenceGenerator : public IReferenceGenerator {
  /// The reference surface onto which to project
  std::shared_ptr<Surface> referenceSurface = nullptr;

  /// @brief  The number of segments for the polyhedron approximation
  int nSegements = 1;

  /// Absolute expansion value for the reference points
  double expansionValue = 0.0;

  /// Luminous region sampling points for the projection - beam spot
  std::vector<Vector3> luminousRegion = {Vector3(0., 0., -200.),
                                         Vector3(0., 0., 200.0)};

  /// Helper to access the Center point of for filling the grid
  ///
  /// @param gctx the geometry context of this operation
  /// @param surface the surface for which the reference point is to be accessed
  ///
  /// @return a vector of reference points for filling
  const std::vector<Vector3> references(const GeometryContext& gctx,
                                        const Surface& surface) const override;
};

}  // namespace Acts

// This file is part of the Acts project.
//
// Copyright (C) 2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// BoundaryIntersectionSorter.h, Acts project
///////////////////////////////////////////////////////////////////

#pragma once

#include <vector>
#include "Acts/Geometry/BoundarySurfaceT.hpp"
#include "Acts/Geometry/TrackingVolume.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/Intersection.hpp"

namespace Acts {

// Full intersection with surface
using BoundaryIntersection =
    FullIntersection<BoundarySurfaceT<TrackingVolume>, Surface>;

// Typedef of the surface intersection
using SurfaceIntersection = ObjectIntersection<Surface>;

/// @brief This struct sorts the boundary surfaces of a tracking volume. The
/// sorting is based on the intersection with a straight line along propagation
/// direction.
struct DefaultBoundaryIntersectionSorter {
  /// @brief Default constructor
  DefaultBoundaryIntersectionSorter() = default;

  /// @brief Main call operator. It tries to intersect every boundary surface by
  /// a straight line along the propagation direction and returns them sorted by
  /// intersecion probability.
  ///
  /// @tparam options_t Type of navigation options object for decomposition
  /// @tparam corrector_t Type of (optional) corrector for surface intersection
  ///
  /// @param gctx The current geometry context object, e.g. alignment
  /// @param boundaries Vector of boundaries of the volume
  /// @param position The position for searching
  /// @param direction The direction for searching
  /// @param options The templated navigation options
  /// @param corrfnc is the corrector struct / function
  ///
  /// @return Vector of intersections with the boundaries ordered by the
  /// intersection probability
  template <typename options_t, typename corrector_t>
  std::vector<BoundaryIntersection> operator()(
      const GeometryContext& gctx,
      std::vector<const BoundarySurfaceT<TrackingVolume>*>& boundaries,
      const Vector3D& position, const Vector3D& direction,
      const options_t& options, const corrector_t& corrfnc) const {
    std::vector<BoundaryIntersection> bIntersections;
    for (auto& bSurface : boundaries) {
      const auto& bSurfaceRep = bSurface->surfaceRepresentation();
      // intersect the surface
      SurfaceIntersection bsIntersection =
          bSurfaceRep.surfaceIntersectionEstimate(gctx, position, direction,
                                                  options, corrfnc);
      // check if the intersection is valid, but exlude the on-surface case
      // when requested -- move to intersectionestimate
      if (bsIntersection) {
        bIntersections.push_back(
            BoundaryIntersection(bsIntersection.intersection, bSurface,
                                 &bSurfaceRep, options.navDir));
      }
    }
    // and now sort to get the closest - need custom sort here to respect sign
    // sort them accordingly to the path length
    if (options.navDir == forward) {
      std::sort(bIntersections.begin(), bIntersections.end());
    } else {
      std::sort(bIntersections.begin(), bIntersections.end(), std::greater<>());
    }
    return bIntersections;
  }
};

/// @brief This struct sorts the boundary surfaces of a tracking volume. The
/// sorting is based on the probability of intersection from a given location in
/// the phase space. The ordering itself is split into three part:
/// 1) Order all surfaces in the given direction based on their path length
/// (most probable)
/// 2) Add all non-intersecting surfaces (somewhat probable in e.g. scattering
/// or curvature cases)
/// 3) Order all surfaces in the opposite direction based on their path length
/// (least probable)
struct BoundaryIntersectionSorter {
  /// @brief Constructor
  BoundaryIntersectionSorter() = default;

  /// @brief Main call operator. The function takes all boundaries and orders
  /// them according to their intersection probability (= path length). Beside
  /// the possible interactions in the given direction the surfaces in the
  /// opposite direction are ordered to receive the least probable surfaces.
  /// This provides the possibility to navigate backwards. Each surface without
  /// a straight line intersection in forward or backward direction is assumed
  /// to be (quite) orthogonal to the current direction and therewith assumed to
  /// be somewhat likely to be intersected. Therewith the navigation in magnetic
  /// fields can be ensured.
  ///
  /// @tparam options_t Type of navigation options object for decomposition
  /// @tparam corrector_t Type of (optional) corrector for surface intersection
  ///
  /// @param gctx The current geometry context object, e.g. alignment
  /// @param boundaries Vector of boundaries of the volume
  /// @param position The position for searching
  /// @param direction The direction for searching
  /// @param options The templated navigation options
  /// @param corrfnc is the corrector struct / function
  ///
  /// @return Vector of intersections with the boundaries ordered by the
  /// intersection probability
  template <typename options_t, typename corrector_t>
  std::vector<BoundaryIntersection> operator()(
      const GeometryContext& gctx,
      std::vector<const BoundarySurfaceT<TrackingVolume>*>& boundaries,
      const Vector3D& position, const Vector3D& direction,
      const options_t& options, const corrector_t& corrfnc) const {
    // Resulting vector
    std::vector<BoundaryIntersection> bIntersections, bIntersectionsOtherNavDir;
    bIntersections.reserve(boundaries.size());

    // Set options for forward and backward direction
    options_t optionsOtherNavDir(options);
    switch (options.navDir) {
      case forward: {
        optionsOtherNavDir.navDir = backward;
        break;
      }
      case backward: {
        optionsOtherNavDir.navDir = forward;
        break;
      }
      // anyDirection does not have any opposite direction so it will directly
      // collected
      case anyDirection: {
        // Collect intersections
        for (const auto& boundary : boundaries) {
          bIntersections.push_back(
              intersect(gctx, boundary, position, direction, options, corrfnc));
        }
        // Fast exit for that case
        return bIntersections;
      }
    }

    // Collect the boundary surfaces both directions
    for (const auto& boundary : boundaries) {
      // If the intersection is valid it is pushed to the final vector otherwise
      // to a tmp storage
      BoundaryIntersection intersection =
          intersect(gctx, boundary, position, direction, options, corrfnc);
      if (intersection) {
        bIntersections.push_back(std::move(intersection));
      } else {
        bIntersectionsOtherNavDir.push_back(intersect(
            gctx, boundary, position, direction, optionsOtherNavDir, corrfnc));
      }
    }

    // Sort both lists
    if (options.navDir == forward) {
      std::sort(bIntersections.begin(), bIntersections.end());
      std::sort(bIntersectionsOtherNavDir.begin(),
                bIntersectionsOtherNavDir.end(), std::greater<>());
    } else {
      std::sort(bIntersections.begin(), bIntersections.end(), std::greater<>());
      std::sort(bIntersectionsOtherNavDir.begin(),
                bIntersectionsOtherNavDir.end());
    }

    // @p bIntersectionsOtherNavDir is sorted such that the backwards elements
    // are in front and then the non-interacting elements. So, they get inserted
    // in inverted order to preserve the probability ordering in the original
    // direction.
    bIntersections.insert(bIntersections.end(),
                          bIntersectionsOtherNavDir.rbegin(),
                          bIntersectionsOtherNavDir.rend());

    return bIntersections;
  }

 private:
  /// @brief This function constructs the intersection with a boundary surface
  /// and returns the intersection object
  ///
  /// @tparam options_t Type of navigation options object for decomposition
  /// @tparam corrector_t Type of (optional) corrector for surface intersection
  ///
  /// @param gctx The current geometry context object, e.g. alignment
  /// @param boundary Boundary of the volume
  /// @param position The position for searching
  /// @param direction The direction for searching
  /// @param options The templated navigation options
  /// @param corrfnc is the corrector struct / function
  ///
  /// @return Intersection object of the boundary surface
  template <typename options_t, typename corrector_t>
  BoundaryIntersection intersect(
      const GeometryContext& gctx,
      const BoundarySurfaceT<TrackingVolume>* boundary,
      const Vector3D& position, const Vector3D& direction,
      const options_t& options, const corrector_t& corrfnc) const {
    const Surface* surface = &(boundary->surfaceRepresentation());
    // intersect the surface
    SurfaceIntersection bsIntersection = surface->surfaceIntersectionEstimate(
        gctx, position, direction, options, corrfnc);
    return BoundaryIntersection(bsIntersection.intersection, boundary, surface,
                                options.navDir);
  }
};
}  // namespace Acts

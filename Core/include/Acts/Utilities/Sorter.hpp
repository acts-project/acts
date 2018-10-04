// This file is part of the Acts project.
//
// Copyright (C) 2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// Sorter.h, Acts project
///////////////////////////////////////////////////////////////////

#pragma once

#include <vector>
#include "Acts/Detector/TrackingVolume.hpp"

namespace Acts {

/// @brief This struct sorts the boundary surfaces of a tracking volume. The
/// sorting is based on the probability of intersection from a given location in
/// the phase space. The ordering itself is splitted into three part:
/// 1) Order all surfaces in the given direction based on their path length
/// (most probable)
/// 2) Add all non-intersecting surfaces (somewhat probable in e.g. scattering
/// or curvature cases)
/// 3) Order all surfaces in the opposite direction based on their path length
/// (least probable)
struct BoundaryIntersectionSorter
{
  /// @brief Constructor
  BoundaryIntersectionSorter() = default;

  /// @brief Main call operator. The function takes all boundaries and orders
  /// them according to their intersection probability (= path length). Beside
  /// the possible interactions in the given direction the surfaces in the
  /// opposite direction are ordered to receive the least probable surfaces.
  /// This provides the possibility to navigate backwards. Every surface which
  /// does not have an interaction is assumed to be (quite) orthogonal to the
  /// current direction and therewith assumed to be somewhat likely to be
  /// intersected. Therewith the navigation in magnetic fields can be ensured.
  ///
  /// @tparam parameters_t Type of parameters used for the decomposition
  /// @tparam options_t Type of navigation options object for decomposition
  /// @tparam corrector_t Type of (optional) corrector for surface intersection
  ///
  /// @param boundaries Vector of boundaries of the volume
  /// @param parameters The templated parameters for searching
  /// @param options The templated navigation options
  /// @param corrfnc is the corrector struct / function
  ///
  /// @return Vector of intersections with the boundaries ordered by the
  /// intersection probability
  template <typename parameters_t, typename options_t, typename corrector_t>
  std::vector<BoundaryIntersection>
  operator()(std::vector<const BoundarySurfaceT<TrackingVolume>*>& boundaries,
             const parameters_t&                                   parameters,
             const options_t&                                      options,
             const corrector_t& corrfnc) const
  {
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
      std::vector<const BoundarySurfaceT<TrackingVolume>*>::iterator
          boundaryIter
          = boundaries.begin();
      for (; boundaryIter != boundaries.end(); boundaryIter++) {
        bIntersections.push_back(
            std::move(intersect(*boundaryIter, parameters, options, corrfnc)));
      }
      // Fast exit for that case
      return bIntersections;
    }
    }

    // Collect the boundary surfaces both directions
    std::vector<const BoundarySurfaceT<TrackingVolume>*>::iterator boundaryIter
        = boundaries.begin();
    for (; boundaryIter != boundaries.end(); boundaryIter++) {
      // If the intersection is valid it is pushed to the final vector otherwise
      // to a tmp storage
      BoundaryIntersection intersection
          = intersect(*boundaryIter, parameters, options, corrfnc);
      if (intersection) {
        bIntersections.push_back(std::move(intersection));
      } else {
        bIntersectionsOtherNavDir.push_back(std::move(
            intersect(*boundaryIter, parameters, optionsOtherNavDir, corrfnc)));
      }
    }

    // Sort both lists
    if (options.navDir == forward) {
      std::sort(bIntersections.begin(), bIntersections.end());
      std::sort(bIntersectionsOtherNavDir.begin(),
                bIntersectionsOtherNavDir.end(),
                std::greater<>());
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
  /// @tparam parameters_t Type of parameters used for the decomposition
  /// @tparam options_t Type of navigation options object for decomposition
  /// @tparam corrector_t Type of (optional) corrector for surface intersection
  ///
  /// @param boundary Boundary of the volume
  /// @param parameters The templated parameters for searching
  /// @param options The templated navigation options
  /// @param corrfnc is the corrector struct / function
  ///
  /// @return Intersection object of the boundary surface
  template <typename parameters_t, typename options_t, typename corrector_t>
  BoundaryIntersection
  intersect(const BoundarySurfaceT<TrackingVolume>* boundary,
            const parameters_t&                     parameters,
            const options_t&                        options,
            const corrector_t&                      corrfnc) const
  {
    Surface const* surface = &(boundary->surfaceRepresentation());
    // intersect the surface
    SurfaceIntersection bsIntersection
        = surface->intersectionEstimate(parameters, options, corrfnc);
    return BoundaryIntersection(
        bsIntersection.intersection, boundary, surface, options.navDir);
  }
};

}  // end of Acts namespace

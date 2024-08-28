// This file is part of the Acts project.
//
// Copyright (C) 2024 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Detector/Portal.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Navigation/NavigationStream.hpp"
#include "Acts/Surfaces/BoundaryTolerance.hpp"
#include "Acts/Utilities/Intersection.hpp"

#include <vector>

namespace Acts {

using namespace Experimental;

class Surface;

namespace NavigationStreamHelper {

/// Helper struct that allows to fill a surface into the candidate vector
///
/// @param nStream the navigation stream that is being filled
/// @param surface the surface to be filled
/// @param bTolerance the boundary tolerance used for the intersection
/// @param abortTarget a boolean that indicates if this surface is an abort target
inline static void fillSurface(NavigationStream& nStream,
                               const Surface* surface,
                               BoundaryTolerance bTolerance,
                               bool abortTarget = false) {
  nStream.candidates.push_back(
      NavigationStream::Candidate{ObjectIntersection<Surface>::invalid(surface),
                                  nullptr, bTolerance, abortTarget});
}

/// Helper struct that allows to fill surfaces into the candidate vector
///
/// @param nStream the navigation stream that is being filled
/// @param surfaces the surfaces that are filled in
/// @param bTolerance the boundary tolerance used for the intersection
/// @param abortTarget a boolean that indicates if those surfaces are abort targets
inline static void fillSurfaces(NavigationStream& nStream,
                                const std::vector<const Surface*>& surfaces,
                                BoundaryTolerance bTolerance,
                                bool abortTarget = false) {
  std::for_each(surfaces.begin(), surfaces.end(), [&](const auto& s) {
    nStream.candidates.push_back(
        NavigationStream::Candidate{ObjectIntersection<Surface>::invalid(s),
                                    nullptr, bTolerance, abortTarget});
  });
}

/// Helper struct that allows to fill portals into the candidate vector
///
/// @param nStream the navigation stream that is being filled
/// @param portals the portals that are filled in
inline static void fillPortals(NavigationStream& nStream,
                               const std::vector<const Portal*>& portals) {
  std::for_each(portals.begin(), portals.end(), [&](const auto& p) {
    nStream.candidates.push_back(NavigationStream::Candidate{
        ObjectIntersection<Surface>::invalid(&(p->surface())), p,
        BoundaryTolerance::None()});
  });
}

/// Initialize a stream that does not require a state object
///
/// @param stream [in, out] is the navigation stream to be updated
/// @param gctx is the geometry context
/// @param queryPoint holds current position, direction, etc.
/// @param cTolerance is the candidate search tolerance
/// @param onSurfaceTolerance is the tolerance for on-surface intersections
///
/// This method will first de-duplicate the candidates on basis of the surface
/// pointer to make sure that the multi-intersections are handled correctly.
/// This will allow intializeStream() to be called even as a re-initialization
/// and still work correctly with at one time valid candidates.
///
/// @return true if the stream is active, false indicates that there are no valid candidates
bool initializeStream(NavigationStream& stream, const GeometryContext& gctx,
                      const NavigationStream::QueryPoint& queryPoint,
                      BoundaryTolerance cTolerance,
                      ActsScalar onSurfaceTolerance = s_onSurfaceTolerance);

/// Convenience method to update a stream from a new query point,
/// this could be called from navigation delegates that do not require
/// a local state or from the navigator on the target stream
///
/// @param stream [in, out] is the navigation stream to be updated
/// @param gctx is the geometry context
/// @param queryPoint holds current position, direction, etc.
/// @param onSurfaceTolerance is the tolerance for on-surface intersections
///
/// @return true if the stream is active, false indicate no valid candidates left
bool updateStream(NavigationStream& stream, const GeometryContext& gctx,
                  const NavigationStream::QueryPoint& queryPoint,
                  ActsScalar onSurfaceTolerance = s_onSurfaceTolerance);

}  // namespace NavigationStreamHelper

}  // namespace Acts

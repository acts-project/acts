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

/// Fillers and attachers for surfaces to act on the navigation state
class NavigationStreamHelper {
 public:
  /// Helper struct that allows to fill surfaces into the candidate vector it
  /// allows to use common navigation structs for volume, portal, surfaces
  ///
  /// @param nStream the navigation stream that is being filled
  /// @param surfaces the surfaces that are filled in
  /// @param bTolerance the boundary tolerance used for the intersection
  inline static void fillSurfaces(NavigationStream& nStream,
                                  const std::vector<const Surface*>& surfaces,
                                  BoundaryTolerance bTolerance) {
    std::for_each(surfaces.begin(), surfaces.end(), [&](const auto& s) {
      nStream.candidates.push_back(NavigationStream::Candidate{
          ObjectIntersection<Surface>(s), nullptr, bTolerance});
    });
  }

  /// Helper struct that allows to fill surfaces into the candidate vector it
  /// allows to use common navigation structs for volume, portal, surfaces
  ///
  /// @param nStream the navigation stream that is being filled
  /// @param portals the portals that are filled in
  inline static void fillPortals(NavigationStream& nStream,
                                 const std::vector<const Portal*>& portals) {
    std::for_each(portals.begin(), portals.end(), [&](const auto& p) {
      nStream.candidates.push_back(NavigationStream::Candidate{
          ObjectIntersection<Surface>(&(p->surface())), p,
          BoundaryTolerance::None()});
    });
  }

  /// Initialize a stream that does not require a state object
  ///
  /// @param stream is the navigation stream to be updated
  /// @param gctx is the geometry context
  /// @param position is the current position
  /// @param direction is the current direction
  /// @param cTolerance is the candidate search tolerance
  ///
  /// @return true if the stream is active, false indicates that there are no valid candidates
  inline static bool initializeStream(NavigationStream& stream,
                                      const GeometryContext& gctx,
                                      const Vector3& position,
                                      const Vector3& direction,
                                      BoundaryTolerance cTolerance) {
    return processStream(stream, gctx, position, direction, true, cTolerance);
  }

  /// Convenience method to upsate a stream from a new position,
  /// this could be called from navigation delegates that do not require
  /// a local state or from the navigator on the target stream
  ///
  /// @param stream is the navigation stream to be updated
  /// @param gctx is the geometry context
  /// @param position is the current position
  /// @param direction is the current direction
  ///
  /// @return true if the stream is active, false indicate no valid candidates left
  inline static bool updateStream(NavigationStream& stream,
                                  const GeometryContext& gctx,
                                  const Vector3& position,
                                  const Vector3& direction) {
    return processStream(stream, gctx, position, direction, false,
                         BoundaryTolerance::None());
  }

 private:
  /// Private helper method that serves intialization and update for
  ///
  /// @param stream is the navigation stream to be updated
  /// @param gctx is the geometry context
  /// @param position is the current position
  /// @param direction is the current direction
  /// @param initial is the flag for the initial update
  /// @param cTolerance is the candidate search tolerance (used only for initial updates)
  static bool processStream(NavigationStream& stream,
                            const GeometryContext& gctx,
                            const Vector3& position, const Vector3& direction,
                            bool initial, BoundaryTolerance cTolerance);

};  // namespace NavigationStreamHelper

}  // namespace Acts

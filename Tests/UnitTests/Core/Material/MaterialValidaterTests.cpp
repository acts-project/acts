// This file is part of the Acts project.
//
// Copyright (C) 2024 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "Acts/Material/HomogeneousSurfaceMaterial.hpp"
#include "Acts/Material/MaterialInteraction.hpp"
#include "Acts/Material/MaterialMapper.hpp"
#include "Acts/Material/MaterialValidater.hpp"
#include "Acts/Material/MaterialSlab.hpp"
#include "Acts/Material/interface/IAssignmentFinder.hpp"
#include "Acts/Surfaces/CylinderSurface.hpp"

#include <limits>

namespace Acts::Test {

auto tContext = GeometryContext();

/// @brief Interface for the material mapping that seeks the possible
/// assignment candidates for the material interactiosn
class IntersectSurfacesFinder : public IAssignmentFinder {
 public:
  std::vector<const Surface*> surfaces;

  /// @brief Interface method for generating assignment candidates for the
  /// material interaction assignment to surfaces or volumes
  ///
  /// @param gctx is the geometry context
  /// @param mctx is the magnetic field context
  /// @param position is the position of the initial ray
  /// @param direction is the direction of initial ray
  ///
  /// @return a vector of Surface Assignments and Volume Assignments
  std::pair<std::vector<IAssignmentFinder::SurfaceAssignment>,
            std::vector<IAssignmentFinder::VolumeAssignment>>
  assignmentCandidates(const GeometryContext& gctx,
                       const MagneticFieldContext& /*ignored*/,
                       const Vector3& position,
                       const Vector3& direction) const override {
    std::vector<IAssignmentFinder::SurfaceAssignment> surfaceAssignments;
    std::vector<IAssignmentFinder::VolumeAssignment> volumeAssignments;
    // Intersect the surfaces
    for (auto& surface : surfaces) {
      // Get the intersection
      auto sMultiIntersection = surface->intersect(gctx, position, direction,
                                                   Acts::BoundaryCheck(true));
      // One solution, take it
      if (sMultiIntersection.size() == 1u &&
          sMultiIntersection[0u].status() >=
              Acts::IntersectionStatus::reachable &&
          sMultiIntersection[0u].pathLength() >= 0.0) {
        surfaceAssignments.push_back(
            {surface, sMultiIntersection[0u].position(), direction});
        continue;
      }
      if (sMultiIntersection.size() > 1u) {
        // Multiple intersections, take the closest
        auto closestForward = sMultiIntersection.closestForward();
        if (closestForward.status() >= Acts::IntersectionStatus::reachable &&
            closestForward.pathLength() > 0.0) {
          surfaceAssignments.push_back(
              {surface, closestForward.position(), direction});
          continue;
        }
      }
    }
    return {surfaceAssignments, volumeAssignments};
  }
};


BOOST_AUTO_TEST_SUITE(MAterialValidatorTestSuite)

BOOST_AUTO_TEST_CASE(MaterialValidaterFlowTest) {

  auto cylinder0 = Surface::makeShared<CylinderSurface>(
      CylinderSurface::Bounds(0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
      CylinderSurface::Cylinder(0.0, 0.0, 0.0, 0.0, 0.0, 0.0));


}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace Acts::Test

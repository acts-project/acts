// This file is part of the Acts project.
//
// Copyright (C) 2017-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/test/data/test_case.hpp>
#include <boost/test/tools/output_test_stream.hpp>
#include <boost/test/unit_test.hpp>

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/Tolerance.hpp"
#include "Acts/Definitions/Units.hpp"
#include "Acts/Geometry/GenericApproachDescriptor.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/Layer.hpp"
#include "Acts/Surfaces/BoundaryCheck.hpp"
#include "Acts/Surfaces/CylinderSurface.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Surfaces/SurfaceArray.hpp"
#include "Acts/Tests/CommonHelpers/FloatComparisons.hpp"
#include "Acts/Utilities/Intersection.hpp"

#include <cstddef>
#include <limits>
#include <memory>
#include <vector>

#include "../Surfaces/SurfaceStub.hpp"
#include "LayerStub.hpp"

using boost::test_tools::output_test_stream;
namespace utf = boost::unit_test;

namespace Acts {

namespace Test {
namespace Layers {

// Build a default context for testing
GeometryContext tgContext = GeometryContext();

BOOST_AUTO_TEST_SUITE(Layers)

/// Unit test for creating compliant/non-compliant GenericApproachDescriptor
/// object
BOOST_AUTO_TEST_CASE(GenericApproachDescriptorConstruction) {
  std::vector<std::shared_ptr<const Surface>> someSurfaces{
      Surface::makeShared<SurfaceStub>(), Surface::makeShared<SurfaceStub>()};
  BOOST_CHECK_NO_THROW(
      GenericApproachDescriptor minimallyConstructedApproachDescriptor(
          someSurfaces));
  //
  std::vector<std::shared_ptr<const Layer>> sharedLayers{
      std::make_shared<LayerStub>(nullptr),
      std::make_shared<LayerStub>(nullptr)};
  BOOST_CHECK_NO_THROW(GenericApproachDescriptor sharedLayerApproachDescriptor(
      {sharedLayers.at(0)->surfaceRepresentation().getSharedPtr(),
       sharedLayers.at(1)->surfaceRepresentation().getSharedPtr()}));
}

/// Unit test for testing GenericApproachDescriptor properties
BOOST_AUTO_TEST_CASE(GenericApproachDescriptorProperties) {
  Vector3 origin{
      0.,
      0.,
      0.,
  };
  Vector3 zDir{0., 0., 1.};
  BoundaryCheck bcheck{true};
  double pLimit = std::numeric_limits<double>::max();
  double oLimit = -100 * UnitConstants::um;
  double tolerance = s_onSurfaceTolerance;
  //
  std::vector<std::shared_ptr<const Surface>> someSurfaces{
      Surface::makeShared<SurfaceStub>(), Surface::makeShared<SurfaceStub>()};
  GenericApproachDescriptor approachDescriptor(someSurfaces);
  LayerStub aLayer(nullptr);
  // registerLayer()
  BOOST_CHECK_NO_THROW(approachDescriptor.registerLayer(aLayer));
  // approachSurface
  SurfaceIntersection surfIntersection = approachDescriptor.approachSurface(
      tgContext, origin, zDir, bcheck, pLimit, oLimit, tolerance);
  double expectedIntersection = 20.0;  // property of SurfaceStub
  CHECK_CLOSE_REL(surfIntersection.pathLength(), expectedIntersection, 1e-6);
  // containedSurfaces()
  BOOST_CHECK_EQUAL(approachDescriptor.containedSurfaces().size(),
                    someSurfaces.size());

  for (size_t i = 0; i < someSurfaces.size(); i++) {
    BOOST_CHECK_EQUAL(approachDescriptor.containedSurfaces().at(i),
                      someSurfaces.at(i).get());
  }
}

/// Unit test for testing GenericApproachDescriptor overstepping
/// - for the approach estimate, there is no overstepping tolerance
/// allowed
BOOST_AUTO_TEST_CASE(GenericApproachNoOverstepping) {
  Vector3 origin{0., -0.5, 1.};
  Vector3 direction{0., 1., 0.};
  BoundaryCheck bcheck{true};
  double pLimit = std::numeric_limits<double>::max();
  double oLimit = -100 * UnitConstants::um;
  double tolerance = s_onSurfaceTolerance;

  auto conCyl =
      Surface::makeShared<CylinderSurface>(Transform3::Identity(), 10., 20.);

  std::vector<std::shared_ptr<const Surface>> approachSurface = {conCyl};

  GenericApproachDescriptor gad(approachSurface);

  auto sfIntersection = gad.approachSurface(
      GeometryContext(), origin, direction, bcheck, pLimit, oLimit, tolerance);

  // No overstepping allowed, the preferred solution should be the forward one
  CHECK_CLOSE_ABS(sfIntersection.pathLength(), 10.5, s_epsilon);
  CHECK_CLOSE_ABS(sfIntersection.position().x(), 0., s_epsilon);
  CHECK_CLOSE_ABS(sfIntersection.position().y(), 10., s_epsilon);
  CHECK_CLOSE_ABS(sfIntersection.position().z(), 1., s_epsilon);
}

BOOST_AUTO_TEST_SUITE_END()
}  // namespace Layers
}  // namespace Test

}  // namespace Acts

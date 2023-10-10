// This file is part of the Acts project.
//
// Copyright (C) 2017-2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/test/data/test_case.hpp>
#include <boost/test/tools/output_test_stream.hpp>
#include <boost/test/unit_test.hpp>

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Surfaces/LineBounds.hpp"
#include "Acts/Surfaces/RectangleBounds.hpp"
#include "Acts/Surfaces/StrawSurface.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Tests/CommonHelpers/DetectorElementStub.hpp"

#include <memory>
#include <string>

namespace Acts {
class PlanarBounds;
}  // namespace Acts

using boost::test_tools::output_test_stream;
namespace utf = boost::unit_test;

namespace Acts {

namespace Test {

// Create a test context
GeometryContext tgContext = GeometryContext();

BOOST_AUTO_TEST_SUITE(StrawSurfaces)
/// Unit test for creating compliant/non-compliant StrawSurface object
BOOST_AUTO_TEST_CASE(StrawSurfaceConstruction) {
  // StrawSurface default constructor is deleted
  //
  /// Constructor with transform, radius and halfZ
  double radius(1.0), halfZ(10.);
  Translation3 translation{0., 1., 2.};
  auto pTransform = Transform3(translation);
  BOOST_CHECK_EQUAL(
      Surface::makeShared<StrawSurface>(Transform3::Identity(), radius, halfZ)
          ->type(),
      Surface::Straw);
  BOOST_CHECK_EQUAL(
      Surface::makeShared<StrawSurface>(pTransform, radius, halfZ)->type(),
      Surface::Straw);
  //
  /// Constructor with transform and LineBounds pointer
  auto pLineBounds = std::make_shared<const LineBounds>(radius, halfZ);
  BOOST_CHECK_EQUAL(
      Surface::makeShared<StrawSurface>(pTransform, pLineBounds)->type(),
      Surface::Straw);
  //
  /// Constructor with LineBounds ptr, DetectorElement
  std::shared_ptr<const Acts::PlanarBounds> p =
      std::make_shared<const RectangleBounds>(1., 10.);
  DetectorElementStub detElement{pTransform, p, 1.0, nullptr};
  BOOST_CHECK_EQUAL(
      Surface::makeShared<StrawSurface>(pLineBounds, detElement)->type(),
      Surface::Straw);
  //
  /// Copy constructor
  auto strawSurfaceObject =
      Surface::makeShared<StrawSurface>(pTransform, radius, halfZ);
  auto copiedStrawSurface =
      Surface::makeShared<StrawSurface>(*strawSurfaceObject);
  BOOST_CHECK_EQUAL(copiedStrawSurface->type(), Surface::Straw);
  BOOST_CHECK(*copiedStrawSurface == *strawSurfaceObject);
  //
  /// Copied and transformed
  auto copiedTransformedStrawSurface = Surface::makeShared<StrawSurface>(
      tgContext, *strawSurfaceObject, pTransform);
  BOOST_CHECK_EQUAL(copiedTransformedStrawSurface->type(), Surface::Straw);
}
//
/// Unit test for testing StrawSurface properties
BOOST_AUTO_TEST_CASE(StrawSurfaceProperties) {
  /// Test clone method
  double radius(1.0), halfZ(10.);
  Translation3 translation{0., 1., 2.};
  auto pTransform = Transform3(translation);
  auto strawSurfaceObject =
      Surface::makeShared<StrawSurface>(pTransform, radius, halfZ);
  //
  /// Test type (redundant)
  BOOST_CHECK_EQUAL(strawSurfaceObject->type(), Surface::Straw);
  //
  /// Test name
  BOOST_CHECK_EQUAL(strawSurfaceObject->name(),
                    std::string("Acts::StrawSurface"));
  //
  /// Test dump
  boost::test_tools::output_test_stream dumpOuput;
  strawSurfaceObject->toStream(tgContext, dumpOuput);
  BOOST_CHECK(
      dumpOuput.is_equal("Acts::StrawSurface\n\
     Center position  (x, y, z) = (0.0000, 1.0000, 2.0000)\n\
     Rotation:             colX = (1.000000, 0.000000, 0.000000)\n\
                           colY = (0.000000, 1.000000, 0.000000)\n\
                           colZ = (0.000000, 0.000000, 1.000000)\n\
     Bounds  : Acts::LineBounds: (radius, halflengthInZ) = (1.0000000, 10.0000000)"));
}

BOOST_AUTO_TEST_CASE(EqualityOperators) {
  double radius(1.0), halfZ(10.);
  Translation3 translation{0., 1., 2.};
  auto pTransform = Transform3(translation);
  auto strawSurfaceObject =
      Surface::makeShared<StrawSurface>(pTransform, radius, halfZ);
  //
  auto strawSurfaceObject2 =
      Surface::makeShared<StrawSurface>(pTransform, radius, halfZ);
  //
  /// Test equality operator
  BOOST_CHECK(*strawSurfaceObject == *strawSurfaceObject2);
  //
  BOOST_TEST_CHECKPOINT(
      "Create and then assign a StrawSurface object to the existing one");
  /// Test assignment
  auto assignedStrawSurface =
      Surface::makeShared<StrawSurface>(Transform3::Identity(), 6.6, 33.33);
  *assignedStrawSurface = *strawSurfaceObject;
  /// Test equality of assigned to original
  BOOST_CHECK(*assignedStrawSurface == *strawSurfaceObject);
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace Test

}  // namespace Acts

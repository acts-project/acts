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

#include <limits>

#include "Acts/Surfaces/ConeSurface.hpp"
#include "Acts/Surfaces/RectangleBounds.hpp"
#include "Acts/Tests/CommonHelpers/FloatComparisons.hpp"
#include "Acts/Utilities/Definitions.hpp"

namespace tt = boost::test_tools;
using boost::test_tools::output_test_stream;
namespace utf = boost::unit_test;

namespace Acts {

namespace Test {

// Create a test context
GeometryContext tgContext = GeometryContext();

BOOST_AUTO_TEST_SUITE(ConeSurfaces)
/// Unit test for creating compliant/non-compliant ConeSurface object
BOOST_AUTO_TEST_CASE(ConeSurfaceConstruction) {
  // ConeSurface default constructor is deleted
  //
  /// Constructor with transform pointer, null or valid, alpha and symmetry
  /// indicator
  double alpha{M_PI / 8.}, halfPhiSector{M_PI / 16.}, zMin{1.0}, zMax{10.};
  bool symmetric(false);
  Translation3D translation{0., 1., 2.};
  auto pTransform = std::make_shared<const Transform3D>(translation);
  std::shared_ptr<const Transform3D> pNullTransform{};
  BOOST_CHECK_EQUAL(
      Surface::makeShared<ConeSurface>(pNullTransform, alpha, symmetric)
          ->type(),
      Surface::Cone);
  BOOST_CHECK_EQUAL(
      Surface::makeShared<ConeSurface>(pTransform, alpha, symmetric)->type(),
      Surface::Cone);
  //
  /// Constructor with transform pointer, alpha,z min and max, halfPhiSector
  BOOST_CHECK_EQUAL(Surface::makeShared<ConeSurface>(pTransform, alpha, zMin,
                                                     zMax, halfPhiSector)
                        ->type(),
                    Surface::Cone);
  //

  /// Constructor with transform and ConeBounds pointer
  // ConeBounds (double alpha, double zmin, double zmax, double halfphi=M_PI,
  // double avphi=0.)
  auto pConeBounds =
      std::make_shared<const ConeBounds>(alpha, zMin, zMax, halfPhiSector, 0.);
  BOOST_CHECK_EQUAL(
      Surface::makeShared<ConeSurface>(pTransform, pConeBounds)->type(),
      Surface::Cone);
  //
  //
  /// Copy constructor
  auto coneSurfaceObject =
      Surface::makeShared<ConeSurface>(pTransform, alpha, symmetric);
  auto copiedConeSurface = Surface::makeShared<ConeSurface>(*coneSurfaceObject);
  BOOST_CHECK_EQUAL(copiedConeSurface->type(), Surface::Cone);
  BOOST_CHECK(*copiedConeSurface == *coneSurfaceObject);
  //
  /// Copied and transformed
  auto copiedTransformedConeSurface = Surface::makeShared<ConeSurface>(
      tgContext, *coneSurfaceObject, *pTransform);
  BOOST_CHECK_EQUAL(copiedTransformedConeSurface->type(), Surface::Cone);

  /// Construct with nullptr bounds
  BOOST_CHECK_THROW(
      auto nullBounds = Surface::makeShared<ConeSurface>(nullptr, nullptr),
      AssertionFailureException);
}
//
/// Unit test for testing ConeSurface properties
BOOST_AUTO_TEST_CASE(ConeSurfaceProperties) {
  /// Test clone method
  double alpha{M_PI / 8.} /*,halfPhiSector{M_PI/16.}, zMin{1.0}, zMax{10.}*/;
  bool symmetric(false);
  Translation3D translation{0., 1., 2.};
  auto pTransform = std::make_shared<const Transform3D>(translation);
  auto coneSurfaceObject =
      Surface::makeShared<ConeSurface>(pTransform, alpha, symmetric);
  //
  /// Test type (redundant)
  BOOST_CHECK_EQUAL(coneSurfaceObject->type(), Surface::Cone);
  //
  /// Test binningPosition
  Vector3D binningPosition{0., 1., 2.};
  CHECK_CLOSE_ABS(
      coneSurfaceObject->binningPosition(tgContext, BinningValue::binPhi),
      binningPosition, 1e-6);
  //
  /// Test referenceFrame
  Vector3D globalPosition{2.0, 2.0, 2.0};
  Vector3D momentum{1.e6, 1.e6, 1.e6};
  double rootHalf = std::sqrt(0.5);
  RotationMatrix3D expectedFrame;
  expectedFrame << -rootHalf, 0., rootHalf, rootHalf, 0., rootHalf, 0., 1., 0.;
  CHECK_CLOSE_OR_SMALL(
      coneSurfaceObject->referenceFrame(tgContext, globalPosition, momentum),
      expectedFrame, 1e-6, 1e-9);
  //
  /// Test normal, given 3D position
  Vector3D origin{0., 0., 0.};
  Vector3D normal3D = {0., -1., 0.};
  CHECK_CLOSE_ABS(coneSurfaceObject->normal(tgContext, origin), normal3D, 1e-6);
  //
  /// Test normal given 2D rphi position
  Vector2D positionPiBy2(1.0, M_PI / 2.);
  Vector3D normalAtPiBy2{0.0312768, 0.92335, -0.382683};

  CHECK_CLOSE_OR_SMALL(coneSurfaceObject->normal(tgContext, positionPiBy2),
                       normalAtPiBy2, 1e-2, 1e-9);
  //
  /// Test rotational symmetry axis
  Vector3D symmetryAxis{0., 0., 1.};
  CHECK_CLOSE_ABS(coneSurfaceObject->rotSymmetryAxis(tgContext), symmetryAxis,
                  1e-6);
  //
  /// Test bounds
  BOOST_CHECK_EQUAL(coneSurfaceObject->bounds().type(), SurfaceBounds::eCone);
  //
  /// Test localToGlobal
  Vector2D localPosition{1.0, M_PI / 2.0};
  coneSurfaceObject->localToGlobal(tgContext, localPosition, momentum,
                                   globalPosition);
  // std::cout<<globalPosition<<std::endl;
  Vector3D expectedPosition{0.0220268, 1.65027, 3.5708};

  CHECK_CLOSE_REL(globalPosition, expectedPosition, 1e-2);
  //
  /// Testing globalToLocal
  coneSurfaceObject->globalToLocal(tgContext, globalPosition, momentum,
                                   localPosition);
  // std::cout<<localPosition<<std::endl;
  Vector2D expectedLocalPosition{1.0, M_PI / 2.0};

  CHECK_CLOSE_REL(localPosition, expectedLocalPosition, 1e-6);
  //
  /// Test isOnSurface
  Vector3D offSurface{100, 1, 2};
  BOOST_CHECK(coneSurfaceObject->isOnSurface(tgContext, globalPosition,
                                             momentum, true));
  BOOST_CHECK(
      !coneSurfaceObject->isOnSurface(tgContext, offSurface, momentum, true));

  /// Test pathCorrection
  CHECK_CLOSE_REL(coneSurfaceObject->pathCorrection(tgContext, offSurface,
                                                    momentum.normalized()),
                  0.40218866453252877, 0.01);
  //
  /// Test name
  BOOST_CHECK_EQUAL(coneSurfaceObject->name(),
                    std::string("Acts::ConeSurface"));
  //
  /// Test dump
  // TODO 2017-04-12 msmk: check how to correctly check output
  //    boost::test_tools::output_test_stream dumpOuput;
  //    coneSurfaceObject.toStream(dumpOuput);
  //    BOOST_CHECK(dumpOuput.is_equal(
  //      "Acts::ConeSurface\n"
  //      "    Center position  (x, y, z) = (0.0000, 1.0000, 2.0000)\n"
  //      "    Rotation:             colX = (1.000000, 0.000000, 0.000000)\n"
  //      "                          colY = (0.000000, 1.000000, 0.000000)\n"
  //      "                          colZ = (0.000000, 0.000000, 1.000000)\n"
  //      "    Bounds  : Acts::ConeBounds: (tanAlpha, minZ, maxZ, averagePhi,
  //      halfPhiSector) = (0.4142136, 0.0000000, inf, 0.0000000,
  //      3.1415927)"));
}

BOOST_AUTO_TEST_CASE(EqualityOperators) {
  double alpha{M_PI / 8.} /*, halfPhiSector{M_PI/16.}, zMin{1.0}, zMax{10.}*/;
  bool symmetric(false);
  Translation3D translation{0., 1., 2.};
  auto pTransform = std::make_shared<const Transform3D>(translation);
  auto coneSurfaceObject =
      Surface::makeShared<ConeSurface>(pTransform, alpha, symmetric);
  //
  auto coneSurfaceObject2 =
      Surface::makeShared<ConeSurface>(pTransform, alpha, symmetric);
  //
  /// Test equality operator
  BOOST_CHECK(*coneSurfaceObject == *coneSurfaceObject2);
  //
  BOOST_TEST_CHECKPOINT(
      "Create and then assign a ConeSurface object to the existing one");
  /// Test assignment
  auto assignedConeSurface =
      Surface::makeShared<ConeSurface>(nullptr, 0.1, true);
  *assignedConeSurface = *coneSurfaceObject;
  /// Test equality of assigned to original
  BOOST_CHECK(*assignedConeSurface == *coneSurfaceObject);
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace Test

}  // namespace Acts

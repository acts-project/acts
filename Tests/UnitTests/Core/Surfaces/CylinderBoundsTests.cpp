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

#include "Acts/Surfaces/CylinderBounds.hpp"
#include "Acts/Tests/CommonHelpers/FloatComparisons.hpp"
#include "Acts/Utilities/Definitions.hpp"

namespace Acts {

namespace Test {
BOOST_AUTO_TEST_SUITE(Surfaces)
/// Unit test for creating compliant/non-compliant CylinderBounds object
BOOST_AUTO_TEST_CASE(CylinderBoundsConstruction) {
  /// test default construction
  // CylinderBounds defaultConstructedCylinderBounds;  // deleted
  double radius(0.5), halfz(10.), halfphi(M_PI / 2.0), averagePhi(M_PI / 2.0);
  BOOST_CHECK_EQUAL(CylinderBounds(radius, halfz).type(),
                    SurfaceBounds::eCylinder);
  BOOST_CHECK_EQUAL(CylinderBounds(radius, halfphi, halfz).type(),
                    SurfaceBounds::eCylinder);
  BOOST_CHECK_EQUAL(CylinderBounds(radius, averagePhi, halfphi, halfz).type(),
                    SurfaceBounds::eCylinder);
  //
  /// test copy construction;
  CylinderBounds cylinderBounds(radius, halfz);
  CylinderBounds copyConstructedCylinderBounds(cylinderBounds);
  BOOST_CHECK_EQUAL(copyConstructedCylinderBounds.type(),
                    SurfaceBounds::eCylinder);
}

/// Unit tests for CylinderBounds properties
BOOST_AUTO_TEST_CASE_EXPECTED_FAILURES(CylinderBoundsProperties, 4)
BOOST_AUTO_TEST_CASE(CylinderBoundsProperties) {
  // CylinderBounds object of radius 0.5 and halfz 20
  double nominalRadius{0.5};
  double nominalHalfLength{20.};
  double averagePhi(0.0);
  double halfphi(M_PI / 4.0);
  CylinderBounds cylinderBoundsObject(nominalRadius, nominalHalfLength);
  CylinderBounds cylinderBoundsSegment(nominalRadius, averagePhi, halfphi,
                                       nominalHalfLength);
  /// test for clone
  auto pCylinderBoundsClone = cylinderBoundsObject.clone();
  BOOST_CHECK_NE(pCylinderBoundsClone, nullptr);
  delete pCylinderBoundsClone;

  /// test for type()
  BOOST_CHECK_EQUAL(cylinderBoundsObject.type(), SurfaceBounds::eCylinder);

  /// test for inside(), 2D coords are r or phi ,z? : needs clarification
  const Vector2D origin{0., 0.};
  const Vector2D atPiBy2{M_PI / 2., 0.0};
  const Vector2D atPi{M_PI, 0.0};
  const Vector2D beyondEnd{0, 30.0};
  const Vector2D unitZ{0.0, 1.0};
  const Vector2D unitPhi{1.0, 0.0};
  const BoundaryCheck trueBoundaryCheckWithTolerance(true, true, 0.1, 0.1);
  BOOST_CHECK(
      cylinderBoundsObject.inside(atPiBy2, trueBoundaryCheckWithTolerance));
  BOOST_CHECK(
      !cylinderBoundsSegment.inside(unitPhi, trueBoundaryCheckWithTolerance));
  BOOST_CHECK(
      cylinderBoundsObject.inside(origin, trueBoundaryCheckWithTolerance));

  /// test for inside3D() with Vector3D argument
  const Vector3D origin3D{0., 0., 0.};
  BOOST_CHECK(
      !cylinderBoundsObject.inside3D(origin3D, trueBoundaryCheckWithTolerance));

  /// test for distanceToBoundary
  CHECK_CLOSE_REL(cylinderBoundsObject.distanceToBoundary(origin), 0.5,
                  1e-6);  // fail
  CHECK_CLOSE_REL(cylinderBoundsObject.distanceToBoundary(beyondEnd), 10.0,
                  1e-6);  // pass
  double sinPiBy8 = std::sin(M_PI / 8.);
  CHECK_CLOSE_REL(cylinderBoundsSegment.distanceToBoundary(atPi), sinPiBy8,
                  1e-6);  // pass
  CHECK_CLOSE_REL(cylinderBoundsSegment.distanceToBoundary(origin), 0.5,
                  1e-6);  // fail

  /// test for r()
  CHECK_CLOSE_REL(cylinderBoundsObject.get(CylinderBounds::eRadius),
                  nominalRadius, 1e-6);

  /// test for averagePhi
  CHECK_CLOSE_REL(cylinderBoundsObject.get(CylinderBounds::eAveragePhi),
                  averagePhi, 1e-6);

  /// test for halfPhiSector
  CHECK_CLOSE_REL(cylinderBoundsSegment.get(CylinderBounds::eHalfPhiSector),
                  halfphi,
                  1e-6);  // fail

  /// test for halflengthZ (NOTE: Naming violation)
  CHECK_CLOSE_REL(cylinderBoundsObject.get(CylinderBounds::eHalfLengthZ),
                  nominalHalfLength, 1e-6);

  /// test for dump
  boost::test_tools::output_test_stream dumpOuput;
  cylinderBoundsObject.toStream(dumpOuput);
  BOOST_CHECK(
      dumpOuput.is_equal("Acts::CylinderBounds: (radius, averagePhi, "
                         "halfPhiSector, halflengthInZ) = (0.5000000, "
                         "0.0000000, 3.1415927, 20.0000000)"));
}
/// Unit test for testing CylinderBounds assignment
BOOST_AUTO_TEST_CASE(CylinderBoundsAssignment) {
  double nominalRadius{0.5};
  double nominalHalfLength{20.};
  CylinderBounds cylinderBoundsObject(nominalRadius, nominalHalfLength);
  CylinderBounds assignedCylinderBounds(10.5, 6.6);
  assignedCylinderBounds = cylinderBoundsObject;
  BOOST_CHECK_EQUAL(assignedCylinderBounds.get(CylinderBounds::eRadius),
                    cylinderBoundsObject.get(CylinderBounds::eRadius));
  BOOST_CHECK_EQUAL(assignedCylinderBounds, cylinderBoundsObject);
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace Test

}  // namespace Acts

// This file is part of the ACTS project.
//
// Copyright (C) 2016 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#define BOOST_TEST_MODULE Cylinder Bounds Tests

#include <boost/test/included/unit_test.hpp>
// leave blank line

#include <boost/test/data/test_case.hpp>
// leave blank line

#include <boost/test/output_test_stream.hpp>
// leave blank line

//
#include "ACTS/Surfaces/CylinderBounds.hpp"
#include "ACTS/Utilities/Definitions.hpp"
//
#include <limits>

// namespace bdata = boost::unit_test::data;
namespace utf    = boost::unit_test;
const double inf = std::numeric_limits<double>::infinity();
const double NaN = std::numeric_limits<double>::quiet_NaN();

namespace Acts {

namespace Test {
  BOOST_AUTO_TEST_SUITE(Surfaces);
  /// Unit test for creating compliant/non-compliant CylinderBounds object
  BOOST_AUTO_TEST_CASE(CylinderBoundsConstruction)
  {
    /// test default construction
    // CylinderBounds defaultConstructedCylinderBounds;  // deleted
    double radius(0.5), halfz(10.), halfphi(M_PI / 2.0), averagePhi(M_PI / 2.0);
    BOOST_TEST(CylinderBounds(radius, halfz).type() == SurfaceBounds::Cylinder);
    BOOST_TEST(CylinderBounds(radius, halfphi, halfz).type()
               == SurfaceBounds::Cylinder);
    BOOST_TEST(CylinderBounds(radius, averagePhi, halfphi, halfz).type()
               == SurfaceBounds::Cylinder);
    //
    /// test copy construction;
    CylinderBounds cylinderBounds(radius, halfz);
    CylinderBounds copyConstructedCylinderBounds(cylinderBounds);
    BOOST_TEST(copyConstructedCylinderBounds.type() == SurfaceBounds::Cylinder,
               "Test trivial method on copy constructed object");
  }
  
  ///Unit tests for CylinderBounds properties
  BOOST_AUTO_TEST_CASE_EXPECTED_FAILURES(CylinderBoundsProperties, 4);
  BOOST_AUTO_TEST_CASE(CylinderBoundsProperties)
  {
    // CylinderBounds object of radius 0.5 and halfz 20
    double nominalRadius{0.5};
    // double r=nominalRadius;
    double         nominalHalfLength{20.};
    double         averagePhi(0.0);
    double         halfphi(M_PI / 4.0);
    double         oneTenthTolerance(0.1);
    CylinderBounds cylinderBoundsObject(nominalRadius, nominalHalfLength);
    CylinderBounds cylinderBoundsSegment(
        nominalRadius, averagePhi, halfphi, nominalHalfLength);
    /// test for clone
    auto pCylinderBoundsClone = cylinderBoundsObject.clone();
    BOOST_CHECK(pCylinderBoundsClone != nullptr);
    delete pCylinderBoundsClone;

    /// test for type()
    BOOST_TEST(cylinderBoundsObject.type() == SurfaceBounds::Cylinder);

    /// test for inside(), 2D coords are r or phi ,z? : needs clarification
    const Vector2D      origin{0., 0.};
    const Vector2D      atPiBy2{M_PI / 2., 0.0};
    const Vector2D      atPi{M_PI, 0.0};
    const Vector2D      beyondEnd{0, 30.0};
    const Vector2D      unitZ{0.0, 1.0};
    const Vector2D      unitPhi{1.0, 0.0};
    const BoundaryCheck trueBoundaryCheckWithTolerance(true, true, 0.1, 0.1);
    BOOST_TEST(
        cylinderBoundsObject.inside(atPiBy2, trueBoundaryCheckWithTolerance)
        == true);
    BOOST_TEST(
        cylinderBoundsSegment.inside(unitPhi, trueBoundaryCheckWithTolerance)
        == false);

    BOOST_TEST(
        cylinderBoundsObject.inside(origin, trueBoundaryCheckWithTolerance)
        == true);

    /// test for insideLoc0 and insideLoc1
    BOOST_TEST(cylinderBoundsObject.insideLoc0(origin, oneTenthTolerance)
               == true);
    BOOST_TEST(cylinderBoundsObject.insideLoc0(unitZ, oneTenthTolerance)
               == true);
    BOOST_TEST(cylinderBoundsObject.insideLoc1(unitPhi, oneTenthTolerance)
               == true);
    BOOST_TEST(cylinderBoundsObject.insideLoc0(beyondEnd, oneTenthTolerance)
               == false);

    /// test for inside3D() with Vector3D argument
    const Vector3D origin3D{0., 0., 0.};
    BOOST_TEST(
        cylinderBoundsObject.inside3D(origin3D, trueBoundaryCheckWithTolerance)
        == false);

    /// test for distanceToBoundary
    BOOST_TEST(cylinderBoundsObject.distanceToBoundary(origin) == 0.5);  // fail
    BOOST_TEST(cylinderBoundsObject.distanceToBoundary(beyondEnd)
               == 10.0);  // pass
    double sinPiBy8 = std::sin(M_PI / 8.);
    BOOST_TEST(cylinderBoundsSegment.distanceToBoundary(atPi)
               == sinPiBy8);  // pass
    BOOST_TEST(cylinderBoundsSegment.distanceToBoundary(origin)
               == 0.5);  // fail

    /// test for r()
    BOOST_TEST(cylinderBoundsObject.r() == nominalRadius);

    /// test for averagePhi
    BOOST_TEST(cylinderBoundsObject.averagePhi() == averagePhi);

    /// test for halfPhiSector
    BOOST_TEST(cylinderBoundsSegment.halfPhiSector() == halfphi);  // fail

    /// test for halflengthZ (NOTE: Naming violation)
    BOOST_TEST(cylinderBoundsObject.halflengthZ() == nominalHalfLength);

    /// test for dump
    boost::test_tools::output_test_stream dumpOuput;
    cylinderBoundsObject.dump(dumpOuput);
    BOOST_TEST(dumpOuput.is_equal("Acts::CylinderBounds: (radius, averagePhi, "
                                  "halfPhiSector, halflengthInZ) = (0.5000000, "
                                  "0.0000000, 3.1415927, 20.0000000)"));
  }
  /// Unit test for testing CylinderBounds assignment
  BOOST_AUTO_TEST_CASE(CylinderBoundsAssignment)
  {
    double         nominalRadius{0.5};
    double         nominalHalfLength{20.};
    CylinderBounds cylinderBoundsObject(nominalRadius, nominalHalfLength);
    CylinderBounds assignedCylinderBounds(NaN, NaN);
    assignedCylinderBounds = cylinderBoundsObject;
    BOOST_TEST(assignedCylinderBounds.r() == cylinderBoundsObject.r());
    BOOST_TEST(assignedCylinderBounds == cylinderBoundsObject);
  }
  BOOST_AUTO_TEST_SUITE_END();

}  // end of namespace Test

}  // end of namespace Acts

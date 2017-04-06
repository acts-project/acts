// This file is part of the ACTS project.
//
// Copyright (C) 2016 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#define BOOST_TEST_MODULE Radial Bounds Tests

#include <boost/test/included/unit_test.hpp>
// leave blank

#include <boost/test/data/test_case.hpp>
// leave blank

#include <boost/test/output_test_stream.hpp>
// leave blank

//
#include "ACTS/Surfaces/RadialBounds.hpp"
#include "ACTS/Utilities/Definitions.hpp"
//
#include <limits>

const double NaN = std::numeric_limits<double>::quiet_NaN();

namespace Acts {

namespace Test {
  BOOST_AUTO_TEST_SUITE(Surfaces);
  /// Unit tests for RadialBounds constrcuctors
  BOOST_AUTO_TEST_CASE(RadialBoundsConstruction)
  {
    double minRadius(1.0), maxRadius(5.0), halfPhiSector(M_PI / 8.0);
    // test default construction
    // RadialBounds defaultConstructedRadialBounds;  should be deleted
    //
    /// Test construction with radii and default sector
    BOOST_TEST(RadialBounds(minRadius, maxRadius).type()
               == SurfaceBounds::Disc);
    //
    /// Test construction with radii and sector half angle
    BOOST_TEST(RadialBounds(minRadius, maxRadius, halfPhiSector).type()
               == SurfaceBounds::Disc);
    //
    /// Copy constructor
    RadialBounds original(minRadius, maxRadius);
    RadialBounds copied(original);
    BOOST_TEST(copied.type() == SurfaceBounds::Disc);
  }

  /// Unit tests for RadialBounds properties
  BOOST_AUTO_TEST_CASE(RadialBoundsProperties)
  {
    double minRadius(1.0), maxRadius(5.0), halfPhiSector(M_PI / 8.0);
    /// Test clone
    RadialBounds radialBoundsObject(minRadius, maxRadius, halfPhiSector);
    auto         pClonedRadialBounds = radialBoundsObject.clone();
    BOOST_TEST(bool(pClonedRadialBounds));
    delete pClonedRadialBounds;
    //
    /// Test type() (redundant; already used in constructor confirmation)
    BOOST_TEST(radialBoundsObject.type() == SurfaceBounds::Disc);
    //
    /// Test distanceToBoundary
    Vector2D origin(0., 0.);
    Vector2D outside(30., 0.);
    Vector2D inSurface(2., 0.0);
    BOOST_TEST(radialBoundsObject.distanceToBoundary(origin)
               == 1.);  // makes sense
    BOOST_TEST(radialBoundsObject.distanceToBoundary(outside) == 25.);  // ok
    //
    /// Test dump
    boost::test_tools::output_test_stream dumpOuput;
    radialBoundsObject.dump(dumpOuput);
    BOOST_TEST(dumpOuput.is_equal("Acts::RadialBounds:  (innerRadius, "
                                  "outerRadius, hPhiSector) = (1.0000000, "
                                  "5.0000000, 0.0000000, 0.3926991)"));
    //
    /// Test inside
    BOOST_TEST(radialBoundsObject.inside(inSurface, BoundaryCheck(true))
               == true);
    BOOST_TEST(radialBoundsObject.inside(outside, BoundaryCheck(true))
               == false);
    //
    /// Test rMin
    BOOST_TEST(radialBoundsObject.rMin() == minRadius);
    //
    /// Test rMax
    BOOST_TEST(radialBoundsObject.rMax() == maxRadius);
    //
    /// Test averagePhi (should be a redundant method, this is not configurable)
    BOOST_TEST(radialBoundsObject.averagePhi() == 0.0);
    //
    /// Test halfPhiSector
    BOOST_TEST(radialBoundsObject.halfPhiSector() == halfPhiSector);
  }
  /// Unit test for testing RadialBounds assignment
  BOOST_AUTO_TEST_CASE(RadialBoundsAssignment)
  {
    double       minRadius(1.0), maxRadius(5.0), halfPhiSector(M_PI / 8.0);
    RadialBounds radialBoundsObject(minRadius, maxRadius, halfPhiSector);
    // operator == not implemented in this class
    //
    /// Test assignment
    RadialBounds assignedRadialBoundsObject(
        NaN, NaN);  // invalid object, in some sense
    assignedRadialBoundsObject = radialBoundsObject;
    BOOST_TEST(assignedRadialBoundsObject == radialBoundsObject);
  }
  BOOST_AUTO_TEST_SUITE_END();

}  // end of namespace Test

}  // end of namespace Acts

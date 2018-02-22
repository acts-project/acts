// This file is part of the ACTS project.
//
// Copyright (C) 2017-2018 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#define BOOST_TEST_MODULE DiscTrapezoidal Bounds Tests

#include <boost/test/included/unit_test.hpp>
// leave blank

#include <boost/test/data/test_case.hpp>
// leave blank

#include <boost/test/output_test_stream.hpp>
// leave blank

//
#include "ACTS/Surfaces/DiscTrapezoidalBounds.hpp"
#include "ACTS/Utilities/Definitions.hpp"

#include <limits>

const double NaN = std::numeric_limits<double>::quiet_NaN();

namespace Acts {

namespace Test {
  BOOST_AUTO_TEST_SUITE(Surfaces)
  /// Unit tests for DiscTrapezoidalBounds constrcuctors
  BOOST_AUTO_TEST_CASE(DiscTrapezoidalBoundsConstruction)
  {
    double minHalfX(1.0), maxHalfX(5.0), rMin(2.0), rMax(6.0), averagePhi(0.0),
        stereo(0.1);
    // test default construction
    // DiscTrapezoidalBounds defaultConstructedDiscTrapezoidalBounds;  should be
    // deleted
    //
    /// Test construction with dimensions and default stereo
    BOOST_TEST(
        DiscTrapezoidalBounds(minHalfX, maxHalfX, rMin, rMax, averagePhi).type()
        == SurfaceBounds::DiscTrapezoidal);
    //
    /// Test construction with all dimensions
    BOOST_TEST(DiscTrapezoidalBounds(
                   minHalfX, maxHalfX, rMin, rMax, averagePhi, stereo)
                   .type()
               == SurfaceBounds::DiscTrapezoidal);
    //
    /// Copy constructor
    DiscTrapezoidalBounds original(minHalfX, maxHalfX, rMin, rMax, averagePhi);
    DiscTrapezoidalBounds copied(original);
    BOOST_TEST(copied.type() == SurfaceBounds::DiscTrapezoidal);
  }

  /// Unit tests for DiscTrapezoidalBounds properties
  BOOST_AUTO_TEST_CASE(DiscTrapezoidalBoundsProperties)
  {
    double minHalfX(1.0), maxHalfX(5.0), rMin(2.0), rMax(6.0),
        averagePhi(0.0) /*, stereo(0.1)*/;
    /// Test clone
    DiscTrapezoidalBounds discTrapezoidalBoundsObject(
        minHalfX, maxHalfX, rMin, rMax, averagePhi);
    auto pClonedDiscTrapezoidalBounds = discTrapezoidalBoundsObject.clone();
    BOOST_TEST(bool(pClonedDiscTrapezoidalBounds));
    delete pClonedDiscTrapezoidalBounds;
    //
    /// Test type() (redundant; already used in constructor confirmation)
    BOOST_TEST(discTrapezoidalBoundsObject.type()
               == SurfaceBounds::DiscTrapezoidal);
    //
    /// Test distanceToBoundary
    Vector2D origin(0., 0.);
    Vector2D outside(30., 0.);
    Vector2D inSurface(2., 0.0);
    BOOST_TEST(discTrapezoidalBoundsObject.distanceToBoundary(origin) == 2.0);
    BOOST_TEST(discTrapezoidalBoundsObject.distanceToBoundary(outside) == 24.0);
    //
    /// Test dump
    boost::test_tools::output_test_stream dumpOuput;
    discTrapezoidalBoundsObject.dump(dumpOuput);
    BOOST_TEST(dumpOuput.is_equal(
        "Acts::DiscTrapezoidalBounds:  (innerRadius, outerRadius, hMinX, "
        "hMaxX, hlengthY, hPhiSector, averagePhi, rCenter, stereo) = "
        "(2.0000000, 6.0000000, 1.0000000, 5.0000000, 0.7922870, 0.9851108, "
        "0.0000000, 2.5243378, 0.0000000)"));
    //
    /// Test inside
    BOOST_TEST(
        discTrapezoidalBoundsObject.inside(inSurface, BoundaryCheck(true))
        == true);
    BOOST_TEST(discTrapezoidalBoundsObject.inside(outside, BoundaryCheck(true))
               == false);
    //
    /// Test rMin
    BOOST_TEST(discTrapezoidalBoundsObject.rMin() == rMin);
    //
    /// Test rMax
    BOOST_TEST(discTrapezoidalBoundsObject.rMax() == rMax);
    //
    /// Test averagePhi
    BOOST_TEST(discTrapezoidalBoundsObject.averagePhi() == averagePhi);
    //
    /// Test rCenter (redundant; not configurable)
    BOOST_TEST(discTrapezoidalBoundsObject.rCenter() == 2.5243377989621383);
    //
    /// Test halfPhiSector (redundant; not configurable)
    BOOST_TEST(discTrapezoidalBoundsObject.stereo() == 0.0);
    //
    /// Test minHalflengthX
    BOOST_TEST(discTrapezoidalBoundsObject.minHalflengthX() == minHalfX);
    //
    /// Test maxHalflengthX
    BOOST_TEST(discTrapezoidalBoundsObject.maxHalflengthX() == maxHalfX);
    //
    /// Test halflengthY
    BOOST_TEST(discTrapezoidalBoundsObject.halflengthY()
               == 0.79228699139326131);
  }
  /// Unit test for testing DiscTrapezoidalBounds assignment
  BOOST_AUTO_TEST_CASE(DiscTrapezoidalBoundsAssignment)
  {
    double minHalfX(1.0), maxHalfX(5.0), rMin(2.0), rMax(6.0), averagePhi(0.0),
        stereo(0.1);
    DiscTrapezoidalBounds discTrapezoidalBoundsObject(
        minHalfX, maxHalfX, rMin, rMax, averagePhi, stereo);
    // operator == not implemented in this class
    //
    /// Test assignment
    DiscTrapezoidalBounds assignedDiscTrapezoidalBoundsObject(
        NaN, NaN, NaN, NaN, NaN);  // invalid object, in some sense
    assignedDiscTrapezoidalBoundsObject = discTrapezoidalBoundsObject;
    BOOST_TEST(assignedDiscTrapezoidalBoundsObject
               == discTrapezoidalBoundsObject);
  }
  BOOST_AUTO_TEST_SUITE_END()

}  // end of namespace Test

}  // end of namespace Acts

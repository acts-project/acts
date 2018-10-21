// This file is part of the Acts project.
//
// Copyright (C) 2017-2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

// clang-format off
#define BOOST_TEST_MODULE Radial Bounds Tests
#include <boost/test/included/unit_test.hpp>
#include <boost/test/data/test_case.hpp>
#include <boost/test/output_test_stream.hpp>
// clang-format on

#include <limits>

#include "Acts/Surfaces/RadialBounds.hpp"
#include "Acts/Tests/CommonHelpers/FloatComparisons.hpp"
#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Utilities/VariantData.hpp"

namespace Acts {

namespace Test {
  BOOST_AUTO_TEST_SUITE(Surfaces)
  /// Unit tests for RadialBounds constrcuctors
  BOOST_AUTO_TEST_CASE(RadialBoundsConstruction)
  {
    double minRadius(1.0), maxRadius(5.0), halfPhiSector(M_PI / 8.0);
    // test default construction
    // RadialBounds defaultConstructedRadialBounds;  should be deleted
    //
    /// Test construction with radii and default sector
    BOOST_CHECK_EQUAL(RadialBounds(minRadius, maxRadius).type(),
                      SurfaceBounds::Disc);
    //
    /// Test construction with radii and sector half angle
    BOOST_CHECK_EQUAL(RadialBounds(minRadius, maxRadius, halfPhiSector).type(),
                      SurfaceBounds::Disc);
    //
    /// Copy constructor
    RadialBounds original(minRadius, maxRadius);
    RadialBounds copied(original);
    BOOST_CHECK_EQUAL(copied.type(), SurfaceBounds::Disc);
  }

  /// Unit tests for RadialBounds properties
  BOOST_AUTO_TEST_CASE(RadialBoundsProperties)
  {
    double minRadius(1.0), maxRadius(5.0), halfPhiSector(M_PI / 8.0);
    /// Test clone
    RadialBounds radialBoundsObject(minRadius, maxRadius, halfPhiSector);
    auto         pClonedRadialBounds = radialBoundsObject.clone();
    BOOST_CHECK_NE(pClonedRadialBounds, nullptr);
    delete pClonedRadialBounds;
    //
    /// Test type() (redundant; already used in constructor confirmation)
    BOOST_CHECK_EQUAL(radialBoundsObject.type(), SurfaceBounds::Disc);
    //
    /// Test distanceToBoundary
    Vector2D origin(0., 0.);
    Vector2D outside(30., 0.);
    Vector2D inSurface(2., 0.0);
    CHECK_CLOSE_REL(radialBoundsObject.distanceToBoundary(origin),
                    1.,
                    1e-6);  // makes sense
    CHECK_CLOSE_REL(radialBoundsObject.distanceToBoundary(outside),
                    25.,
                    1e-6);  // ok
    //
    /// Test dump
    boost::test_tools::output_test_stream dumpOuput;
    radialBoundsObject.dump(dumpOuput);
    BOOST_CHECK(dumpOuput.is_equal("Acts::RadialBounds:  (innerRadius, "
                                   "outerRadius, hPhiSector) = (1.0000000, "
                                   "5.0000000, 0.0000000, 0.3926991)"));
    //
    /// Test inside
    BOOST_CHECK(radialBoundsObject.inside(inSurface, BoundaryCheck(true)));
    BOOST_CHECK(!radialBoundsObject.inside(outside, BoundaryCheck(true)));
    //
    /// Test rMin
    BOOST_CHECK_EQUAL(radialBoundsObject.rMin(), minRadius);
    //
    /// Test rMax
    BOOST_CHECK_EQUAL(radialBoundsObject.rMax(), maxRadius);
    //
    /// Test averagePhi (should be a redundant method, this is not configurable)
    BOOST_CHECK_EQUAL(radialBoundsObject.averagePhi(), 0.0);
    //
    /// Test halfPhiSector
    BOOST_CHECK_EQUAL(radialBoundsObject.halfPhiSector(), halfPhiSector);
  }
  /// Unit test for testing RadialBounds assignment
  BOOST_AUTO_TEST_CASE(RadialBoundsAssignment)
  {
    double       minRadius(1.0), maxRadius(5.0), halfPhiSector(M_PI / 8.0);
    RadialBounds radialBoundsObject(minRadius, maxRadius, halfPhiSector);
    // operator == not implemented in this class
    //
    /// Test assignment
    RadialBounds assignedRadialBoundsObject(10.1, 123.);
    assignedRadialBoundsObject = radialBoundsObject;
    BOOST_CHECK_EQUAL(assignedRadialBoundsObject, radialBoundsObject);
  }

  BOOST_AUTO_TEST_CASE(RadialBounds_toVariantData)
  {
    double       rMin = 1, rMax = 5, avgPhi = M_PI / 3., halfPhi = M_PI;
    RadialBounds rad(rMin, rMax, avgPhi, halfPhi);

    variant_data var_rad = rad.toVariantData();
    std::cout << var_rad << std::endl;

    variant_map var_rad_map = boost::get<variant_map>(var_rad);
    BOOST_CHECK_EQUAL(var_rad_map.get<std::string>("type"), "RadialBounds");
    variant_map pl = var_rad_map.get<variant_map>("payload");
    BOOST_CHECK_EQUAL(pl.get<double>("rMin"), rMin);
    BOOST_CHECK_EQUAL(pl.get<double>("rMax"), rMax);
    BOOST_CHECK_EQUAL(pl.get<double>("avgPhi"), avgPhi);
    BOOST_CHECK_EQUAL(pl.get<double>("halfPhi"), halfPhi);

    RadialBounds rad2(var_rad);

    BOOST_CHECK_EQUAL(rad.rMin(), rad2.rMin());
    BOOST_CHECK_EQUAL(rad.rMax(), rad2.rMax());
    BOOST_CHECK_EQUAL(rad.averagePhi(), rad2.averagePhi());
    BOOST_CHECK_EQUAL(rad.halfPhiSector(), rad2.halfPhiSector());
  }

  BOOST_AUTO_TEST_SUITE_END()

}  // namespace Test

}  // namespace Acts

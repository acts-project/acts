// This file is part of the ACTS project.
//
// Copyright (C) 2016 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#define BOOST_TEST_MODULE BoundaryCheck Tests

#include <boost/test/included/unit_test.hpp>
//leave blank line

#include <boost/test/data/test_case.hpp>
//leave blank line

#include <boost/test/output_test_stream.hpp>
//leave blank line


#include "ACTS/Surfaces/BoundaryCheck.hpp"
#include "ACTS/Utilities/Definitions.hpp"
#include "ACTS/Utilities/Units.hpp"
#include <limits>

namespace bdata = boost::unit_test::data;
namespace utf   = boost::unit_test;
const double inf = std::numeric_limits<double>::infinity();
const double NaN = std::numeric_limits<double>::quiet_NaN();

namespace Acts {

namespace Test {
  
  bool
  approximatelyEqual(const double & a, const double & b, const double & tol){
    const double dif0=std::abs(a - b);
    return  (dif0<tol);
  }
  bool equalKdop(const KDOP & kdop0, const KDOP & kdop1){
    return ( kdop0.min == kdop1.min) and (kdop0.max == kdop1.max);
  }
  BOOST_AUTO_TEST_SUITE(Surfaces);
  // See: https://en.wikipedia.org/wiki/Bounding_volume
  //
  /// Unit test for creating compliant/non-compliant BoundaryCheck object
  BOOST_AUTO_TEST_CASE(BoundaryCheckConstruction)
  { 
    //BoundaryCheck boundaryCheckDefault; no default constructor, doesn't compile
    //BOOST_TEST(boundaryCheckDefault.checkLoc1 == false);
    //all constructors simply set member variables
    //so first we check the simplest case
    BoundaryCheck boundaryCheckTrue{true};
    BOOST_TEST(boundaryCheckTrue.checkLoc1 == true);
    //..and then the most complicated
    ActsSymMatrixD<2> cov;
    cov << 0.04, 0, 0, 0.1;
    const double nSigma=1.2;
    BoundaryCheck boundaryCheckWithCovariance{cov, nSigma, true, true};
    BOOST_TEST(boundaryCheckWithCovariance.checkLoc1 == true);
    BoundaryCheck copyConstructedBoundaryCheck(boundaryCheckWithCovariance);
    const auto & originalCovariance = *(boundaryCheckWithCovariance.lCovariance);
    const auto & copiedCovariance = *(copyConstructedBoundaryCheck.lCovariance);
    BOOST_TEST(originalCovariance == copiedCovariance);
    //corner cases (NaN, inf, in tolerance and covariance) are not tested.
  }
  ///Unit test for computing bounding box along arbitrary axes (but here chosen to be x,y)
  BOOST_AUTO_TEST_CASE_EXPECTED_FAILURES( BoundaryCheckComputeKdop, 1 );
  BOOST_AUTO_TEST_CASE(BoundaryCheckComputeKdop){
    //NOTE: There is no check that the axes are normalised, so this will fail
    const Vector2D longXAxis(2.,0.);
    const Vector2D longYAxis(0.,2.);
    std::vector<Vector2D> longAxes{longXAxis, longYAxis};
    //
    const Vector2D xAxis(1.,0.);
    const Vector2D yAxis(0.,1.);
    std::vector<Vector2D> axes{xAxis, yAxis};
    //
    std::vector<Vector2D> convexPolygon{ {3.,2.}, {5.,2.}, {6.,3.}, {3.,4.}, {2.,3.} };
    //NOTE: KDOP is float precision, everything else is double
    //NOTE: There is no check this has the same dimension as the number of axes
    std::vector<KDOP> boundingExtentsAlongAxes{ {0.,0.}, {0.,0.} };
    BoundaryCheck b(true);
    //NOTE: Naming violation. 
    //NOTE: doesn't need to be class member
    b.ComputeKDOP(convexPolygon, axes, boundingExtentsAlongAxes);
    //expected result
    std::vector<KDOP> expectedResult{ {2.,6.}, {2.,4.} };
    //
    bool boundingBoxMeetsExpectations = equalKdop(expectedResult[0],boundingExtentsAlongAxes[0]);
    boundingBoxMeetsExpectations &= equalKdop(expectedResult[1],boundingExtentsAlongAxes[1]);
    BOOST_TEST(boundingBoxMeetsExpectations,"Bounding box correctly calculated for unit axes");
    //
    b.ComputeKDOP(convexPolygon, longAxes, boundingExtentsAlongAxes);
    boundingBoxMeetsExpectations = equalKdop(expectedResult[0],boundingExtentsAlongAxes[0]);
    boundingBoxMeetsExpectations &= equalKdop(expectedResult[1],boundingExtentsAlongAxes[1]);
    BOOST_TEST(boundingBoxMeetsExpectations,"Bounding box correctly calculated for non unit axes");
  }
  ///Unit test for bounding boxes overlap method
  BOOST_AUTO_TEST_CASE(CheckOverlappingBoundsBoxes){
    std::vector<KDOP> box{ {-2.,2}, {-5.,4.}};
    std::vector<KDOP> overlappingBounds{ {-3., -1.}, {-1.,1.} };
    std::vector<KDOP> nonOverlappingBounds{ {2.5,5.}, {-4.,4.}};
    BoundaryCheck b(true);
    //NOTE: Function doesn't need to be a member
    //NOTE: Naming violation
    BOOST_TEST(b.TestKDOPKDOP(box,overlappingBounds));
    BOOST_TEST(not b.TestKDOPKDOP(box,nonOverlappingBounds));
  }
  
  /// Unit test for testing BoundaryCheck FastSinCos accuracy 
  /// with a data set of 20 steps from -pi to pi
  BOOST_DATA_TEST_CASE(BoundaryCheckFastSinCosAccuracy,
                       bdata::xrange( -M_PI, M_PI, 2.*M_PI/20.0)){
    BoundaryCheck b(true);
    //NOTE: Naming violation
    sincosCache sc = b.FastSinCos(sample);
    const double sinValue=std::sin(sample);
    const double cosValue=std::cos(sample);
    const double fastSin=sc.sinC;
    const double fastCos=sc.cosC;
    //NOTE: Use CHECK_SMALL on the difference to check within absolute tolerance
    //BOOST_TEST(a==b, tolerance) checks the relative tolerance is within range
    BOOST_CHECK_SMALL(fastSin - sinValue, 0.001); //claimed max difference in code
    BOOST_CHECK_SMALL(fastCos - cosValue, 0.001); //claimed max difference in code
  }
  ///Unit test for FastArcTan
  BOOST_DATA_TEST_CASE(BoundaryCheckFastArcTanAccuracy,
                       bdata::xrange( 0., 200., 10.)){
    BoundaryCheck b(true);
    //NOTE: Naming violation
    double fastArctan = b.FastArcTan(sample);
    const double atanValue=std::atan(sample);
    //NOTE: Use CHECK_SMALL on the difference to check within absolute tolerance
    //because BOOST_TEST(a==b, tolerance) checks the relative tolerance.
    BOOST_CHECK_SMALL(fastArctan - atanValue, 0.0015); //claimed max difference in code
  }
  
  ///Unit test for conversion of smooth ellipse to 4 + 4*n point polygon
  BOOST_AUTO_TEST_CASE(BoundaryCheckEllipseToPoly){
    ActsSymMatrixD<2> cov;
    //numbers chosen for ease!
    cov << 0.16, 0., 0., 0.81;
    const double nSigma=1.5;  //what if this is -ve, zero, inf or NaN?
    const int pointsPerQuadrant = 1; //what if this is -ve, zero, inf or NaN?
    BoundaryCheck bc{cov, nSigma, true, true};
    std::vector<Vector2D> convertedEllipsePoints=bc.EllipseToPoly(pointsPerQuadrant);
    //points at (0.6,0), (0.,1.35) and, by plugging in y=(1.35/0.6)*x, x=sqrt(0.18), y=2.25*x
    const double minor{0.6}, major{1.35}, midx{std::sqrt(0.18)}, midy{2.25*midx};
    const std::vector<Vector2D> reference{ {minor,0.0},{-minor,0.0},{0.0,major},\
    {0.0,-major}, {midx,midy},{midx,-midy},{-midx,-midy},{-midx,midy} };
    unsigned int refIdx(0);
    for (auto & i:convertedEllipsePoints){
      auto & j=reference[refIdx++];
      const bool convertedEllipsePointMatchesReferencePoint = approximatelyEqual(i[0],j[0],0.001)
       and approximatelyEqual(i[1],j[1],0.001); //absolute tolerance, not relative.
      BOOST_TEST(convertedEllipsePointMatchesReferencePoint,"point "<<refIdx);
    }
  }
  
  ///Unit test for assignment
  BOOST_AUTO_TEST_CASE_EXPECTED_FAILURES( BoundaryCheckAssignment, 1 );
  BOOST_AUTO_TEST_CASE(BoundaryCheckAssignment){
    //
    ActsSymMatrixD<2> cov;
    cov << 0.04, 0, 0, 0.1;
    const double nSigma=1.2;
    BoundaryCheck original{cov, nSigma, false, true};
    BoundaryCheck assigned(true);
    assigned=original;
    BOOST_TEST(assigned.checkLoc0 == original.checkLoc0);
    BOOST_TEST(assigned.checkLoc1 == original.checkLoc1);
    BOOST_TEST(assigned.toleranceLoc0 == original.toleranceLoc0);
    BOOST_TEST(assigned.toleranceLoc1 == original.toleranceLoc1);
    BOOST_TEST(assigned.nSigmas == original.nSigmas);
    BOOST_TEST( (bool(original.lCovariance)), "Original object pointer to covariance should not be null");
    //the following will fail
    BOOST_TEST( (bool(assigned.lCovariance)), "Assigned object pointer to covariance should not be null");
    if (original.lCovariance and assigned.lCovariance){
      auto & originalCovariance{*(original.lCovariance)};
      auto & assignedCovariance{*(assigned.lCovariance)};
      BOOST_TEST(originalCovariance == assignedCovariance );
    }
  }
  BOOST_AUTO_TEST_SUITE_END();
}  // end of namespace Test

}  // end of namespace Acts

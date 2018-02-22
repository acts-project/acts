// This file is part of the ACTS project.
//
// Copyright (C) 2017-2018 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#define BOOST_TEST_MODULE Cone Bounds Tests

#include <boost/test/included/unit_test.hpp>
// leave blank line

#include <boost/test/data/test_case.hpp>
// leave blank line

#include <boost/test/output_test_stream.hpp>
// leave blank line

//
#include "ACTS/Surfaces/ConeBounds.hpp"
#include "ACTS/Utilities/Definitions.hpp"
//
#include <limits>

const double NaN = std::numeric_limits<double>::quiet_NaN();
/* Note on nomenclature:
  alpha = cone opening half angle
  z is the axis of symmetry
  zmin, zmax define limits for truncated cone
  phi is clock angle around cone, with x axis corresponding to phi=0
  Cone segments may be defined with the avphi (central position of segment) and
    halfphi (extent in phi of cone segment either side of the avphi)
  Local coords are z, rphi
*/
namespace Acts {

namespace Test {
  BOOST_AUTO_TEST_SUITE(Surfaces)
  /// Unit test for creating compliant/non-compliant ConeBounds object
  BOOST_AUTO_TEST_CASE(ConeBoundsConstruction)
  {
    // test default construction
    // ConeBounds defaultConstructedConeBounds;  // deleted
    double alpha(M_PI / 8.0), zMin(3.), zMax(6.), halfPhi(M_PI / 4.0),
        averagePhi(0.);
    const bool symmetric(false);
    BOOST_TEST_CHECKPOINT("Four parameter constructor (last two at default)");
    ConeBounds defaultConeBounds(alpha, symmetric);
    BOOST_TEST(defaultConeBounds.type() == SurfaceBounds::Cone);
    BOOST_TEST_CHECKPOINT("Four parameter constructor");
    ConeBounds fourParameterConstructed(alpha, symmetric, halfPhi, averagePhi);
    BOOST_TEST(fourParameterConstructed.type() == SurfaceBounds::Cone);
    BOOST_TEST_CHECKPOINT("Five parameter constructor (last two at default)");
    ConeBounds defaulted5ParamConeBounds(alpha, zMin, zMax);
    BOOST_TEST(defaulted5ParamConeBounds.type() == SurfaceBounds::Cone);
    BOOST_TEST_CHECKPOINT("Five parameter constructor)");
    ConeBounds fiveParamConstructedConeBounds(
        alpha, zMin, zMax, halfPhi, averagePhi);
    BOOST_TEST(fiveParamConstructedConeBounds.type() == SurfaceBounds::Cone);
    BOOST_TEST_CHECKPOINT("Copy constructor");
    ConeBounds copyConstructedConeBounds(fiveParamConstructedConeBounds);
    BOOST_TEST(copyConstructedConeBounds == fiveParamConstructedConeBounds);
    auto pClonedConeBounds = fiveParamConstructedConeBounds.clone();
    BOOST_TEST(*pClonedConeBounds == fiveParamConstructedConeBounds);
    delete pClonedConeBounds;
  }
  /// Unit tests for properties of ConeBounds object
  BOOST_AUTO_TEST_CASE(ConeBoundsProperties)
  {
    double alpha(M_PI / 8.0), zMin(3.), zMax(6.), halfPhi(M_PI / 4.0),
        averagePhi(0.);
    // const bool symmetric(false);
    const Vector2D origin(0, 0);
    const Vector2D somewhere(4., 4.);
    ConeBounds     coneBoundsObject(alpha, zMin, zMax, halfPhi, averagePhi);
    //
    /// test for type (redundant)
    BOOST_TEST(coneBoundsObject.type() == SurfaceBounds::Cone);
    //
    /// test for inside
    BOOST_TEST(coneBoundsObject.inside(origin) == false);
    //
    /// test for distanceToBoundary
    // std::cout << coneBoundsObject.distanceToBoundary(origin) << std::endl;
    //
    /// test for r
    BOOST_TEST(coneBoundsObject.r(zMin) == zMin * std::tan(alpha));
    //
    /// test for tanAlpha
    BOOST_TEST(coneBoundsObject.tanAlpha() == std::tan(alpha));
    //
    /// test for sinAlpha
    BOOST_TEST(coneBoundsObject.sinAlpha() == std::sin(alpha));
    //
    /// test for cosAlpha
    BOOST_TEST(coneBoundsObject.cosAlpha() == std::cos(alpha));
    //
    /// test for alpha
    BOOST_TEST(coneBoundsObject.alpha() == alpha);
    //
    /// test for minZ
    BOOST_TEST(coneBoundsObject.minZ() == zMin);
    //
    /// test for maxZ
    BOOST_TEST(coneBoundsObject.maxZ() == zMax);
    //
    /// test for averagePhi
    BOOST_TEST(coneBoundsObject.halfPhiSector() == halfPhi);
    //
    /// test for dump
    boost::test_tools::output_test_stream dumpOuput;
    coneBoundsObject.dump(dumpOuput);
    BOOST_TEST(dumpOuput.is_equal(
        "Acts::ConeBounds: (tanAlpha, minZ, maxZ, averagePhi, halfPhiSector) = "
        "(0.4142136, 3.0000000, 6.0000000, 0.0000000, 0.7853982)"));
  }
  /// Unit test for testing ConeBounds assignment
  BOOST_AUTO_TEST_CASE(ConeBoundsAssignment)
  {
    double alpha(M_PI / 8.0), zMin(3.), zMax(6.), halfPhi(M_PI / 4.0),
        averagePhi(0.);
    // const bool symmetric(false);
    ConeBounds originalConeBounds(alpha, zMin, zMax, halfPhi, averagePhi);
    ConeBounds assignedConeBounds(NaN, zMin, zMax, halfPhi, averagePhi);
    assignedConeBounds = originalConeBounds;
    BOOST_TEST(assignedConeBounds == originalConeBounds);
  }
  BOOST_AUTO_TEST_SUITE_END()

}  // end of namespace Test

}  // end of namespace Acts

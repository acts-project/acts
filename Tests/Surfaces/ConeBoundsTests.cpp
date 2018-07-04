// This file is part of the Acts project.
//
// Copyright (C) 2017-2018 Acts project team
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
#include "Acts/Surfaces/ConeBounds.hpp"
#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Utilities/VariantData.hpp"
//
#include "../Utilities/TestHelper.hpp"
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
    BOOST_CHECK_CLOSE_FRACTION(
        coneBoundsObject.r(zMin), zMin * std::tan(alpha), 1e-6);
    //
    /// test for tanAlpha
    BOOST_CHECK_CLOSE_FRACTION(
        coneBoundsObject.tanAlpha(), std::tan(alpha), 1e-6);
    //
    /// test for sinAlpha
    BOOST_CHECK_CLOSE_FRACTION(
        coneBoundsObject.sinAlpha(), std::sin(alpha), 1e-6);
    //
    /// test for cosAlpha
    BOOST_CHECK_CLOSE_FRACTION(
        coneBoundsObject.cosAlpha(), std::cos(alpha), 1e-6);
    //
    /// test for alpha
    BOOST_CHECK_CLOSE_FRACTION(coneBoundsObject.alpha(), alpha, 1e-6);
    //
    /// test for minZ
    BOOST_CHECK_CLOSE_FRACTION(coneBoundsObject.minZ(), zMin, 1e-6);
    //
    /// test for maxZ
    BOOST_CHECK_CLOSE_FRACTION(coneBoundsObject.maxZ(), zMax, 1e-6);
    //
    /// test for averagePhi
    BOOST_CHECK_CLOSE_FRACTION(coneBoundsObject.halfPhiSector(), halfPhi, 1e-6);
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

  BOOST_AUTO_TEST_CASE(ConeBounds_toVariantData)
  {
    double alpha = M_PI / 2., zMin = 1, zMax = 5, avgPhi = M_PI / 3.,
           halfPhi = M_PI;
    ConeBounds cone(alpha, zMin, zMax, halfPhi, avgPhi);

    variant_data var_cone = cone.toVariantData();
    std::cout << var_cone << std::endl;

    variant_map var_cone_map = boost::get<variant_map>(var_cone);
    BOOST_TEST(var_cone_map.get<std::string>("type") == "ConeBounds");
    variant_map pl = var_cone_map.get<variant_map>("payload");
    BOOST_CHECK_CLOSE_FRACTION(pl.get<double>("alpha"), alpha, 1e-6);
    BOOST_CHECK_CLOSE_FRACTION(pl.get<double>("zMin"), zMin, 1e-6);
    BOOST_CHECK_CLOSE_FRACTION(pl.get<double>("zMax"), zMax, 1e-6);
    BOOST_CHECK_CLOSE_FRACTION(pl.get<double>("avgPhi"), avgPhi, 1e-6);
    BOOST_CHECK_CLOSE_FRACTION(pl.get<double>("halfPhi"), halfPhi, 1e-6);

    ConeBounds cone2(var_cone);

    BOOST_CHECK_CLOSE_FRACTION(cone.alpha(), cone2.alpha(), 1e-6);
    BOOST_CHECK_CLOSE_FRACTION(cone.minZ(), cone2.minZ(), 1e-6);
    BOOST_CHECK_CLOSE_FRACTION(cone.maxZ(), cone2.maxZ(), 1e-6);
    BOOST_CHECK_CLOSE_FRACTION(cone.averagePhi(), cone2.averagePhi(), 1e-6);
    BOOST_CHECK_CLOSE_FRACTION(
        cone.halfPhiSector(), cone2.halfPhiSector(), 1e-6);
  }

  BOOST_AUTO_TEST_SUITE_END()

}  // namespace Test

}  // namespace Acts

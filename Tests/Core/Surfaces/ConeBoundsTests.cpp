// This file is part of the Acts project.
//
// Copyright (C) 2017-2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

// clang-format off
#define BOOST_TEST_MODULE Cone Bounds Tests
#include <boost/test/included/unit_test.hpp>
#include <boost/test/data/test_case.hpp>
#include <boost/test/output_test_stream.hpp>
// clang-format on

#include <limits>

#include "Acts/Surfaces/ConeBounds.hpp"
#include "Acts/Tests/CommonHelpers/FloatComparisons.hpp"
#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Utilities/VariantData.hpp"

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
    BOOST_CHECK_EQUAL(defaultConeBounds.type(), SurfaceBounds::Cone);
    BOOST_TEST_CHECKPOINT("Four parameter constructor");
    ConeBounds fourParameterConstructed(alpha, symmetric, halfPhi, averagePhi);
    BOOST_CHECK_EQUAL(fourParameterConstructed.type(), SurfaceBounds::Cone);
    BOOST_TEST_CHECKPOINT("Five parameter constructor (last two at default)");
    ConeBounds defaulted5ParamConeBounds(alpha, zMin, zMax);
    BOOST_CHECK_EQUAL(defaulted5ParamConeBounds.type(), SurfaceBounds::Cone);
    BOOST_TEST_CHECKPOINT("Five parameter constructor)");
    ConeBounds fiveParamConstructedConeBounds(
        alpha, zMin, zMax, halfPhi, averagePhi);
    BOOST_CHECK_EQUAL(fiveParamConstructedConeBounds.type(),
                      SurfaceBounds::Cone);
    BOOST_TEST_CHECKPOINT("Copy constructor");
    ConeBounds copyConstructedConeBounds(fiveParamConstructedConeBounds);
    BOOST_CHECK_EQUAL(copyConstructedConeBounds,
                      fiveParamConstructedConeBounds);
    auto pClonedConeBounds = fiveParamConstructedConeBounds.clone();
    BOOST_CHECK_EQUAL(*pClonedConeBounds, fiveParamConstructedConeBounds);
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
    BOOST_CHECK_EQUAL(coneBoundsObject.type(), SurfaceBounds::Cone);
    //
    /// test for inside
    BOOST_CHECK(!coneBoundsObject.inside(origin));
    //
    /// test for distanceToBoundary
    // std::cout << coneBoundsObject.distanceToBoundary(origin) << std::endl;
    //
    /// test for r
    CHECK_CLOSE_REL(coneBoundsObject.r(zMin), zMin * std::tan(alpha), 1e-6);
    //
    /// test for tanAlpha
    CHECK_CLOSE_REL(coneBoundsObject.tanAlpha(), std::tan(alpha), 1e-6);
    //
    /// test for sinAlpha
    CHECK_CLOSE_REL(coneBoundsObject.sinAlpha(), std::sin(alpha), 1e-6);
    //
    /// test for cosAlpha
    CHECK_CLOSE_REL(coneBoundsObject.cosAlpha(), std::cos(alpha), 1e-6);
    //
    /// test for alpha
    CHECK_CLOSE_REL(coneBoundsObject.alpha(), alpha, 1e-6);
    //
    /// test for minZ
    CHECK_CLOSE_REL(coneBoundsObject.minZ(), zMin, 1e-6);
    //
    /// test for maxZ
    CHECK_CLOSE_REL(coneBoundsObject.maxZ(), zMax, 1e-6);
    //
    /// test for averagePhi
    CHECK_CLOSE_REL(coneBoundsObject.halfPhiSector(), halfPhi, 1e-6);
    //
    /// test for dump
    boost::test_tools::output_test_stream dumpOuput;
    coneBoundsObject.dump(dumpOuput);
    BOOST_CHECK(dumpOuput.is_equal(
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
    ConeBounds assignedConeBounds(0.1, 2.3, 4.5, 1.2, 2.1);
    assignedConeBounds = originalConeBounds;
    BOOST_CHECK_EQUAL(assignedConeBounds, originalConeBounds);
  }

  BOOST_AUTO_TEST_CASE(ConeBounds_toVariantData)
  {
    double alpha = M_PI / 2., zMin = 1, zMax = 5, avgPhi = M_PI / 3.,
           halfPhi = M_PI;
    ConeBounds cone(alpha, zMin, zMax, halfPhi, avgPhi);

    variant_data var_cone = cone.toVariantData();
    std::cout << var_cone << std::endl;

    variant_map var_cone_map = boost::get<variant_map>(var_cone);
    BOOST_CHECK_EQUAL(var_cone_map.get<std::string>("type"), "ConeBounds");
    variant_map pl = var_cone_map.get<variant_map>("payload");
    BOOST_CHECK_EQUAL(pl.get<double>("alpha"), alpha);
    BOOST_CHECK_EQUAL(pl.get<double>("zMin"), zMin);
    BOOST_CHECK_EQUAL(pl.get<double>("zMax"), zMax);
    BOOST_CHECK_EQUAL(pl.get<double>("avgPhi"), avgPhi);
    BOOST_CHECK_EQUAL(pl.get<double>("halfPhi"), halfPhi);

    ConeBounds cone2(var_cone);

    BOOST_CHECK_EQUAL(cone.alpha(), cone2.alpha());
    BOOST_CHECK_EQUAL(cone.minZ(), cone2.minZ());
    BOOST_CHECK_EQUAL(cone.maxZ(), cone2.maxZ());
    BOOST_CHECK_EQUAL(cone.averagePhi(), cone2.averagePhi());
    BOOST_CHECK_EQUAL(cone.halfPhiSector(), cone2.halfPhiSector());
  }

  BOOST_AUTO_TEST_SUITE_END()

}  // namespace Test

}  // namespace Acts

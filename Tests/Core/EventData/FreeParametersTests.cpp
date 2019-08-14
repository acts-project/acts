// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

// clang-format off
#define BOOST_TEST_MODULE FreeParameters Tests
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
// clang-format on

#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/EventData/NeutralParameters.hpp"
#include "Acts/Tests/CommonHelpers/FloatComparisons.hpp"
#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Geometry/GeometryContext.hpp"

namespace Acts {
namespace Test {

/// @brief Unit test for free parameters
///
BOOST_AUTO_TEST_CASE(free_initialization) {
  Vector3D pos(0., 1., 2.);
  double t = 3.;
  Vector3D dir(4., 5., 6.);
  double qop = 7.;

  FreeVector params;
  params << pos.x(), pos.y(), pos.z(), t, dir.x(), dir.y(), dir.z(), qop;

  std::unique_ptr<FreeSymMatrix> covPtr = nullptr;

  // Test if the object can be created w/o covariance
  FreeParameters fpwoCov(nullptr, params);
  BOOST_CHECK_EQUAL(fpwoCov.covariance(), nullptr);
  CHECK_CLOSE_ABS(fpwoCov.parameters(), params, 1e-6);

  // Test if the object can be create with covariance
  FreeSymMatrix cov;
  cov << 1., 0., 0., 0., 0., 0., 0., 0., 0., 2., 0., 0., 0., 0., 0., 0., 0., 0.,
      3., 0., 0., 0., 0., 0., 0., 0., 0., 4., 0., 0., 0., 0., 0., 0., 0., 0.,
      5., 0., 0., 0., 0., 0., 0., 0., 0., 6., 0., 0., 0., 0., 0., 0., 0., 0.,
      7., 0., 0., 0., 0., 0., 0., 0., 0., 8.;
  covPtr = std::make_unique<FreeSymMatrix>(cov);
  FreeParameters fp(std::move(covPtr), params);
  CHECK_CLOSE_COVARIANCE(*fp.covariance(), cov, 1e-6);
  CHECK_CLOSE_ABS(fp.parameters(), params, 1e-6);

  // Test == comparison
  BOOST_TEST(fp == fp);
  BOOST_TEST(fp != fpwoCov);

  FreeParameters fpCopyConstr(fp);
  BOOST_TEST(fpCopyConstr == fp);

  FreeParameters fpMoveConstr(
      FreeParameters(std::make_unique<FreeSymMatrix>(cov), params));
  BOOST_TEST(fpMoveConstr == fp);

  FreeParameters* fpCopy = fp.clone();
  BOOST_TEST(*fpCopy == fp);

  // Test copy assignment
  FreeParameters fpCopyAssignment = fp;
  BOOST_TEST(fpCopyAssignment == fp);

  // Test move assignment
  FreeParameters fpMoveAssignment =
      FreeParameters(std::make_unique<FreeSymMatrix>(cov), params);
  BOOST_TEST(fpMoveAssignment == fp);

  /// Repeat constructing and assignment with neutral parameters

  // Test if the object can be created w/o covariance
  NeutralFreeParameters nfpwoCov(nullptr, params);
  BOOST_CHECK_EQUAL(nfpwoCov.covariance(), nullptr);
  CHECK_CLOSE_ABS(nfpwoCov.parameters(), params, 1e-6);

  NeutralFreeParameters nfp(std::make_unique<FreeSymMatrix>(cov), params);
  CHECK_CLOSE_COVARIANCE(*nfp.covariance(), cov, 1e-6);
  CHECK_CLOSE_ABS(nfp.parameters(), params, 1e-6);

  NeutralFreeParameters nfpCopyConstr(nfp);
  BOOST_TEST(nfpCopyConstr == nfp);

  NeutralFreeParameters nfpMoveConstr(
      NeutralFreeParameters(std::make_unique<FreeSymMatrix>(cov), params));
  BOOST_TEST(nfpMoveConstr == nfp);

  NeutralFreeParameters* nfpCopy = nfp.clone();
  BOOST_TEST(*nfpCopy == nfp);

  // Test copy assignment
  NeutralFreeParameters nfpCopyAssignment = nfp;
  BOOST_TEST(nfpCopyAssignment == nfp);

  // Test move assignment
  NeutralFreeParameters nfpMoveAssignment =
      NeutralFreeParameters(std::make_unique<FreeSymMatrix>(cov), params);
  BOOST_TEST(nfpMoveAssignment == nfp);

  /// Test getters/setters

  // Test getter of single elements
  CHECK_CLOSE_ABS(fp.get<0>(), pos.x(), 1e-6);
  CHECK_CLOSE_ABS(fp.get<1>(), pos.y(), 1e-6);
  CHECK_CLOSE_ABS(fp.get<2>(), pos.z(), 1e-6);
  CHECK_CLOSE_ABS(fp.get<3>(), t, 1e-6);
  CHECK_CLOSE_ABS(fp.get<4>(), dir.x(), 1e-6);
  CHECK_CLOSE_ABS(fp.get<5>(), dir.y(), 1e-6);
  CHECK_CLOSE_ABS(fp.get<6>(), dir.z(), 1e-6);
  CHECK_CLOSE_ABS(fp.get<7>(), qop, 1e-6);

  // Test getter of uncertainties
  CHECK_CLOSE_ABS(fp.uncertainty<0>(), std::sqrt(cov(0, 0)), 1e-6);
  CHECK_CLOSE_ABS(fp.uncertainty<1>(), std::sqrt(cov(1, 1)), 1e-6);
  CHECK_CLOSE_ABS(fp.uncertainty<2>(), std::sqrt(cov(2, 2)), 1e-6);
  CHECK_CLOSE_ABS(fp.uncertainty<3>(), std::sqrt(cov(3, 3)), 1e-6);
  CHECK_CLOSE_ABS(fp.uncertainty<4>(), std::sqrt(cov(4, 4)), 1e-6);
  CHECK_CLOSE_ABS(fp.uncertainty<5>(), std::sqrt(cov(5, 5)), 1e-6);
  CHECK_CLOSE_ABS(fp.uncertainty<6>(), std::sqrt(cov(6, 6)), 1e-6);
  CHECK_CLOSE_ABS(fp.uncertainty<7>(), std::sqrt(cov(7, 7)), 1e-6);

  // Test getter of parts of the parameters by their meaning
  CHECK_CLOSE_ABS(fp.position(), pos, 1e-6);
  CHECK_CLOSE_ABS(fp.momentum(), dir / qop, 1e-6);
  CHECK_CLOSE_ABS(fp.charge(), +1., 1e-6);
  BOOST_TEST(nfp.charge() == 0.);
  CHECK_CLOSE_ABS(fp.time(), t, 1e-6);

  // Test setters
  GeometryContext dummy;
  fp.set<0>(dummy, 8.);
  fp.set<1>(dummy, 9.);
  fp.set<2>(dummy, 10.);
  fp.set<3>(dummy, 11.);
  fp.set<4>(dummy, 12.);
  fp.set<5>(dummy, 13.);
  fp.set<6>(dummy, 14.);
  fp.set<7>(dummy, 15.);
  CHECK_CLOSE_ABS(fp.get<0>(), 8., 1e-6);
  CHECK_CLOSE_ABS(fp.get<1>(), 9., 1e-6);
  CHECK_CLOSE_ABS(fp.get<2>(), 10., 1e-6);
  CHECK_CLOSE_ABS(fp.get<3>(), 11., 1e-6);
  CHECK_CLOSE_ABS(fp.get<4>(), 12., 1e-6);
  CHECK_CLOSE_ABS(fp.get<5>(), 13., 1e-6);
  CHECK_CLOSE_ABS(fp.get<6>(), 14., 1e-6);
  CHECK_CLOSE_ABS(fp.get<7>(), 15., 1e-6);
}
}  // namespace Test
}  // namespace Acts
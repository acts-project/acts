// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/EventData/NeutralTrackParameters.hpp"
#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Tests/CommonHelpers/FloatComparisons.hpp"
#include "Acts/Utilities/Definitions.hpp"

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

  std::optional<FreeSymMatrix> cov = std::nullopt;

  // Test if the object can be created w/o covariance
  FreeTrackParameters fpwoCov(cov, params);
  BOOST_CHECK(!fpwoCov.covariance().has_value());
  CHECK_CLOSE_ABS(fpwoCov.parameters(), params, 1e-6);

  // Test if the object can be create with covariance
  *cov << 1., 0., 0., 0., 0., 0., 0., 0., 0., 2., 0., 0., 0., 0., 0., 0., 0.,
      0., 3., 0., 0., 0., 0., 0., 0., 0., 0., 4., 0., 0., 0., 0., 0., 0., 0.,
      0., 5., 0., 0., 0., 0., 0., 0., 0., 0., 6., 0., 0., 0., 0., 0., 0., 0.,
      0., 7., 0., 0., 0., 0., 0., 0., 0., 0., 8.;
  std::optional<FreeSymMatrix> covCpy = *cov;

  FreeTrackParameters fp(covCpy, params);
  CHECK_CLOSE_COVARIANCE(*fp.covariance(), *cov, 1e-6);
  CHECK_CLOSE_ABS(fp.parameters(), params, 1e-6);

  // Test == comparison
  BOOST_CHECK_EQUAL(fp, fp);
  BOOST_CHECK_NE(fp, fpwoCov);

  FreeTrackParameters fpCopyConstr(fp);
  BOOST_CHECK_EQUAL(fpCopyConstr, fp);

  covCpy = *cov;
  FreeTrackParameters fpMoveConstr(FreeTrackParameters(covCpy, params));
  BOOST_CHECK_EQUAL(fpMoveConstr, fp);

  // Test copy assignment
  FreeTrackParameters fpCopyAssignment = fp;
  BOOST_CHECK_EQUAL(fpCopyAssignment, fp);

  // Test move assignment
  covCpy = *cov;
  FreeTrackParameters fpMoveAssignment = FreeTrackParameters(covCpy, params);
  BOOST_CHECK_EQUAL(fpMoveAssignment, fp);

  /// Repeat constructing and assignment with neutral parameters

  // Test if the object can be created w/o covariance
  NeutralFreeTrackParameters nfpwoCov(std::nullopt, params);
  BOOST_CHECK(!nfpwoCov.covariance().has_value());
  CHECK_CLOSE_ABS(nfpwoCov.parameters(), params, 1e-6);

  covCpy = *cov;
  NeutralFreeTrackParameters nfp(covCpy, params);
  CHECK_CLOSE_COVARIANCE(*nfp.covariance(), *cov, 1e-6);
  CHECK_CLOSE_ABS(nfp.parameters(), params, 1e-6);

  NeutralFreeTrackParameters nfpCopyConstr(nfp);
  BOOST_CHECK_EQUAL(nfpCopyConstr, nfp);

  covCpy = *cov;
  NeutralFreeTrackParameters nfpMoveConstr(
      NeutralFreeTrackParameters(covCpy, params));
  BOOST_CHECK_EQUAL(nfpMoveConstr, nfp);

  // Test copy assignment
  NeutralFreeTrackParameters nfpCopyAssignment = nfp;
  BOOST_CHECK_EQUAL(nfpCopyAssignment, nfp);

  // Test move assignment
  covCpy = *cov;
  NeutralFreeTrackParameters nfpMoveAssignment =
      NeutralFreeTrackParameters(covCpy, params);
  BOOST_CHECK_EQUAL(nfpMoveAssignment, nfp);

  /// Test getters/setters

  // Test getter of single elements
  CHECK_CLOSE_ABS(fp.get<eFreePos0>(), pos.x(), 1e-6);
  CHECK_CLOSE_ABS(fp.get<eFreePos1>(), pos.y(), 1e-6);
  CHECK_CLOSE_ABS(fp.get<eFreePos2>(), pos.z(), 1e-6);
  CHECK_CLOSE_ABS(fp.get<eFreeTime>(), t, 1e-6);
  CHECK_CLOSE_ABS(fp.get<eFreeDir0>(), dir.x(), 1e-6);
  CHECK_CLOSE_ABS(fp.get<eFreeDir1>(), dir.y(), 1e-6);
  CHECK_CLOSE_ABS(fp.get<eFreeDir2>(), dir.z(), 1e-6);
  CHECK_CLOSE_ABS(fp.get<eFreeQOverP>(), qop, 1e-6);

  // Test getter of uncertainties
  CHECK_CLOSE_ABS(fp.uncertainty<eFreePos0>(), std::sqrt((*cov)(0, 0)), 1e-6);
  CHECK_CLOSE_ABS(fp.uncertainty<eFreePos1>(), std::sqrt((*cov)(1, 1)), 1e-6);
  CHECK_CLOSE_ABS(fp.uncertainty<eFreePos2>(), std::sqrt((*cov)(2, 2)), 1e-6);
  CHECK_CLOSE_ABS(fp.uncertainty<eFreeTime>(), std::sqrt((*cov)(3, 3)), 1e-6);
  CHECK_CLOSE_ABS(fp.uncertainty<eFreeDir0>(), std::sqrt((*cov)(4, 4)), 1e-6);
  CHECK_CLOSE_ABS(fp.uncertainty<eFreeDir1>(), std::sqrt((*cov)(5, 5)), 1e-6);
  CHECK_CLOSE_ABS(fp.uncertainty<eFreeDir2>(), std::sqrt((*cov)(6, 6)), 1e-6);
  CHECK_CLOSE_ABS(fp.uncertainty<eFreeQOverP>(), std::sqrt((*cov)(7, 7)), 1e-6);

  // Test getter of parts of the parameters by their meaning
  CHECK_CLOSE_ABS(fp.position(), pos, 1e-6);
  CHECK_CLOSE_ABS(fp.momentum(), dir / qop, 1e-6);
  CHECK_CLOSE_ABS(fp.charge(), +1., 1e-6);
  BOOST_CHECK_EQUAL(nfp.charge(), 0.);
  CHECK_CLOSE_ABS(fp.time(), t, 1e-6);
}

}  // namespace Test
}  // namespace Acts

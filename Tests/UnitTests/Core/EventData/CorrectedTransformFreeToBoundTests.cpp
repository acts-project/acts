// This file is part of the Acts project.
//
// Copyright (C) 2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/test/data/test_case.hpp>
#include <boost/test/unit_test.hpp>

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/Definitions/Units.hpp"
#include "Acts/EventData/detail/CorrectedTransformationFreeToBound.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Tests/CommonHelpers/FloatComparisons.hpp"
#include "Acts/Utilities/Intersection.hpp"
#include "Acts/Utilities/UnitVectors.hpp"

#include <algorithm>
#include <cmath>
#include <memory>
#include <optional>
#include <ostream>
#include <string>
#include <tuple>

using namespace Acts;
using namespace Acts::UnitLiterals;

namespace {
constexpr ActsScalar eps = 0.01;
}

BOOST_AUTO_TEST_CASE(CorrectedFreeToBoundTrackParameters) {
  GeometryContext geoCtx;

  const auto loc0 = 0.0;
  const auto loc1 = 0.0;
  const auto phi = 0.0;
  const auto theta = M_PI / 4;
  const auto qOverP = 1 / 1_GeV;
  const auto t = 1_ns;

  const auto resLoc0 = 0.0;
  const auto resLoc1 = 0.0;
  const auto resPhi = 0.25;
  const auto resTheta = 0.25;
  const auto resQOverP = 0.01 / 1_GeV;
  const auto resTime = 0.01_ns;

  // construct two parallel plane surfaces with normal in x direction
  ActsScalar distance = 10_mm;
  auto sSurface =
      Surface::makeShared<PlaneSurface>(Vector3(0, 0, 0), Vector3::UnitX());
  auto eSurface = Surface::makeShared<PlaneSurface>(Vector3(distance, 0, 0),
                                                    Vector3::UnitX());

  // the bound parameters at the starting plane
  BoundVector sBoundParams = BoundVector::Zero();
  sBoundParams << loc0, loc1, phi, theta, qOverP, t;

  // the bound parameters covariance at the starting  plane
  BoundSquareMatrix sBoundCov = BoundSquareMatrix::Zero();
  sBoundCov(eBoundLoc0, eBoundLoc0) = resLoc0 * resLoc0;
  sBoundCov(eBoundLoc1, eBoundLoc1) = resLoc1 * resLoc1;
  sBoundCov(eBoundPhi, eBoundPhi) = resPhi * resPhi;
  sBoundCov(eBoundTheta, eBoundTheta) = resTheta * resTheta;
  sBoundCov(eBoundQOverP, eBoundQOverP) = resQOverP * resQOverP;
  sBoundCov(eBoundTime, eBoundTime) = resTime * resTime;

  Vector3 dir = makeDirectionFromPhiTheta(phi, theta);

  // the intersection of the track with the end surface
  SurfaceIntersection intersection =
      eSurface->intersect(geoCtx, Vector3(0, 0, 0), dir, true).closest();
  Vector3 tpos = intersection.position();
  auto s = intersection.pathLength();

  BOOST_CHECK_EQUAL(s, distance * std::sqrt(2));

  // construct the free parameters vector
  FreeVector eFreeParams = FreeVector::Zero();
  eFreeParams.segment<3>(0) = tpos;
  eFreeParams[eFreeTime] = t;
  eFreeParams.segment<3>(4) = dir;
  eFreeParams[eFreeQOverP] = qOverP;

  // the jacobian from local to global at the starting position
  BoundToFreeMatrix boundToFreeJac =
      sSurface->boundToFreeJacobian(geoCtx, sBoundParams);

  // the transport jacobian without B field
  FreeMatrix transportJac = FreeMatrix::Identity();
  transportJac(eFreePos0, eFreeDir0) = s;
  transportJac(eFreePos1, eFreeDir1) = s;
  transportJac(eFreePos2, eFreeDir2) = s;

  // the free covariance at the start position
  FreeSquareMatrix sFreeCov =
      boundToFreeJac * sBoundCov * boundToFreeJac.transpose();
  // the free covariance at the end position
  FreeSquareMatrix eFreeCov =
      transportJac * sFreeCov * transportJac.transpose();

  // convert free parameters to bound parameters with non-linear correction

  BOOST_TEST_INFO("Transform free parameters vector onto surface "
                  << eSurface->name());

  // the corrected transformation
  auto freeToBoundCorrection = FreeToBoundCorrection(true);
  BOOST_CHECK(freeToBoundCorrection);

  auto transformer =
      detail::CorrectedFreeToBoundTransformer(freeToBoundCorrection);
  auto correctedRes = transformer(eFreeParams, eFreeCov, *eSurface, geoCtx);

  BOOST_CHECK(correctedRes.has_value());
  auto correctedValue = correctedRes.value();
  BoundVector eCorrectedBoundParams = std::get<BoundVector>(correctedValue);
  BoundSquareMatrix eCorrectedBoundCov =
      std::get<BoundSquareMatrix>(correctedValue);

  // the loc0, phi are the same as that without correction
  BOOST_CHECK_EQUAL(eCorrectedBoundParams[eBoundLoc0], loc0);
  BOOST_CHECK_EQUAL(eCorrectedBoundParams[eBoundPhi], phi);
  CHECK_CLOSE_REL(eCorrectedBoundParams[eBoundLoc1], 11.2563, eps);

  BOOST_TEST_INFO("Corrected Bound Params: \n" << eCorrectedBoundParams);
  BOOST_TEST_INFO("Corrected Bound Covariance: \n" << eCorrectedBoundCov);
  // the error for loc0 is the same as that without correction:
  // loc0 at end position = distance * tan(theta) * sin(phi),
  // dloc0/dphi = distance * tan(theta) * cos(phi) = distance,
  // resolution of loc0 at end position = dloc0/dphi * resLoc0 = 2.5
  CHECK_CLOSE_REL(eCorrectedBoundCov(eBoundLoc0, eBoundLoc0), pow(2.5, 2), eps);
}

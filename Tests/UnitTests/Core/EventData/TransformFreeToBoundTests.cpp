// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/test/data/test_case.hpp>
#include <boost/test/unit_test.hpp>
#include <cmath>
#include <limits>

#include "Acts/EventData/detail/TransformationFreeToBound.hpp"
#include "Acts/Tests/CommonHelpers/FloatComparisons.hpp"
#include "Acts/Utilities/UnitVectors.hpp"
#include "Acts/Utilities/Units.hpp"
#include "TrackParametersTestData.hpp"

using namespace Acts;
using namespace Acts::UnitLiterals;

namespace {
constexpr auto eps = std::numeric_limits<BoundParametersScalar>::epsilon();
}

BOOST_AUTO_TEST_SUITE(TransformFreeToBound)

BOOST_DATA_TEST_CASE(
    GlobalToBoundParameters,
    surfaces* posSymmetric* posSymmetric* ts* phis* thetas* qOverPs, surface,
    l0, l1, time, phiInput, theta, qOverP) {
  // phi is ill-defined in forward/backward tracks
  auto phi = ((0 < theta) and (theta < M_PI)) ? phiInput : 0.0;

  GeometryContext geoCtx;
  Vector2D loc(l0, l1);
  Vector3D pos = Vector3D::Zero();
  Vector3D dir = makeDirectionUnitFromPhiTheta(phi, theta);
  // transform reference position
  surface->localToGlobal(geoCtx, loc, dir, pos);

  // convert to free parameters
  BoundVector bv = detail::transformFreeToBoundParameters(
      pos, time, dir, qOverP, *surface, geoCtx);

  BOOST_TEST_INFO("Using surface " << surface->name());
  CHECK_CLOSE_OR_SMALL(bv[eBoundLoc0], l0, eps, eps);
  CHECK_CLOSE_OR_SMALL(bv[eBoundLoc1], l1, eps, eps);
  CHECK_CLOSE_OR_SMALL(bv[eBoundTime], time, eps, eps);
  CHECK_CLOSE_OR_SMALL(bv[eBoundPhi], phi, eps, eps);
  CHECK_CLOSE_OR_SMALL(bv[eBoundTheta], theta, eps, eps);
  CHECK_CLOSE_OR_SMALL(bv[eBoundQOverP], qOverP, eps, eps);
}

BOOST_DATA_TEST_CASE(GlobalToCurvilinearParameters, ts* phis* thetas* qOverPs,
                     time, phiInput, theta, qOverP) {
  // phi is ill-defined in forward/backward tracks
  auto phi = ((0 < theta) and (theta < M_PI)) ? phiInput : 0.0;

  GeometryContext geoCtx;
  Vector3D dir = makeDirectionUnitFromPhiTheta(phi, theta);

  // convert to free parameters
  BoundVector bv =
      detail::transformFreeToCurvilinearParameters(time, dir, qOverP);

  CHECK_SMALL(bv[eBoundLoc0], eps);
  CHECK_SMALL(bv[eBoundLoc1], eps);
  CHECK_CLOSE_OR_SMALL(bv[eBoundTime], time, eps, eps);
  CHECK_CLOSE_OR_SMALL(bv[eBoundPhi], phi, eps, eps);
  CHECK_CLOSE_OR_SMALL(bv[eBoundTheta], theta, eps, eps);
  CHECK_CLOSE_OR_SMALL(bv[eBoundQOverP], qOverP, eps, eps);
}

BOOST_AUTO_TEST_SUITE_END()

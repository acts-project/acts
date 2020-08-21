// This file is part of the Acts project.
//
// Copyright (C) 2017-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/EventData/NeutralTrackParameters.hpp"
#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/Tests/CommonHelpers/FloatComparisons.hpp"
#include "Acts/Utilities/UnitVectors.hpp"
#include "Acts/Utilities/detail/periodic.hpp"

#include <limits>

#include "TrackParametersTestData.hpp"

using namespace Acts;

static constexpr auto eps =
    8 * std::numeric_limits<BoundParametersScalar>::epsilon();

BOOST_AUTO_TEST_SUITE(CurvilinearTrackParameters)

BOOST_DATA_TEST_CASE(
    ConstructCharged,
    posSymmetric* posSymmetric* posSymmetric* ts* phis* thetas*(ps ^ qs ^
                                                                qOverPs),
    x, y, z, time, phiInput, theta, p, q, qOverP) {
  // phi is ill-defined in forward/backward tracks
  const auto phi = ((0 < theta) and (theta < M_PI)) ? phiInput : 0.0;

  const GeometryContext geoCtx;
  const Vector3D pos(x, y, z);
  const Vector3D dir = makeDirectionUnitFromPhiTheta(phi, theta);

  CurvilinearParameters params(std::nullopt, pos, p * dir, q, time);

  CHECK_SMALL(params.get<eBoundLoc0>(), eps);
  CHECK_SMALL(params.get<eBoundLoc1>(), eps);
  CHECK_CLOSE_OR_SMALL(params.get<eBoundTime>(), time, eps, eps);
  CHECK_CLOSE_OR_SMALL(detail::radian_sym(params.get<eBoundPhi>()),
                       detail::radian_sym(phi), eps, eps);
  CHECK_CLOSE_OR_SMALL(params.get<eBoundTheta>(), theta, eps, eps);
  CHECK_CLOSE_OR_SMALL(params.get<eBoundQOverP>(), qOverP, eps, eps);

  CHECK_CLOSE_OR_SMALL(params.position(geoCtx), pos, eps, eps);
  CHECK_CLOSE_OR_SMALL(params.position(), pos, eps, eps);
  CHECK_CLOSE_OR_SMALL(params.time(), time, eps, eps);
  CHECK_CLOSE_OR_SMALL(params.momentum(), p * dir, eps, eps);
  CHECK_CLOSE_REL(params.charge(), q, eps);

  CHECK_CLOSE_OR_SMALL(params.referenceSurface().center(geoCtx), pos, eps, eps);
  CHECK_CLOSE_OR_SMALL(params.referenceSurface().normal(geoCtx), dir, eps, eps);

  // TODO verify reference frame
}

BOOST_DATA_TEST_CASE(
    ConstructNeutral,
    posSymmetric* posSymmetric* posSymmetric* ts* phis* thetas* ps, x, y, z,
    time, phiInput, theta, p) {
  // phi is ill-defined in forward/backward tracks
  const auto phi = ((0 < theta) and (theta < M_PI)) ? phiInput : 0.0;

  const GeometryContext geoCtx;
  const Vector3D pos(x, y, z);
  const Vector3D dir = makeDirectionUnitFromPhiTheta(phi, theta);

  NeutralCurvilinearTrackParameters params(std::nullopt, pos, p * dir, time);

  CHECK_SMALL(params.get<eBoundLoc0>(), eps);
  CHECK_SMALL(params.get<eBoundLoc1>(), eps);
  CHECK_CLOSE_OR_SMALL(params.get<eBoundTime>(), time, eps, eps);
  CHECK_CLOSE_OR_SMALL(detail::radian_sym(params.get<eBoundPhi>()),
                       detail::radian_sym(phi), eps, eps);
  CHECK_CLOSE_OR_SMALL(params.get<eBoundTheta>(), theta, eps, eps);
  CHECK_CLOSE_OR_SMALL(params.get<eBoundQOverP>(), 1 / p, eps, eps);

  CHECK_CLOSE_OR_SMALL(params.position(geoCtx), pos, eps, eps);
  CHECK_CLOSE_OR_SMALL(params.position(), pos, eps, eps);
  CHECK_CLOSE_OR_SMALL(params.time(), time, eps, eps);
  CHECK_CLOSE_OR_SMALL(params.momentum(), p * dir, eps, eps);
  CHECK_SMALL(params.charge(), eps);

  CHECK_CLOSE_OR_SMALL(params.referenceSurface().center(geoCtx), pos, eps, eps);
  CHECK_CLOSE_OR_SMALL(params.referenceSurface().normal(geoCtx), dir, eps, eps);

  // TODO verify reference frame
}

BOOST_AUTO_TEST_SUITE_END()

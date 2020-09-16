// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/test/data/test_case.hpp>
#include <boost/test/unit_test.hpp>

#include "Acts/EventData/detail/TransformationFreeToBound.hpp"
#include "Acts/Tests/CommonHelpers/FloatComparisons.hpp"
#include "Acts/Utilities/UnitVectors.hpp"
#include "Acts/Utilities/Units.hpp"

#include <cmath>
#include <limits>

#include "TrackParametersDatasets.hpp"

using namespace Acts;
using namespace Acts::UnitLiterals;

namespace {
constexpr auto eps = std::numeric_limits<BoundScalar>::epsilon();
}

BOOST_AUTO_TEST_SUITE(TransformFreeToBound)

BOOST_DATA_TEST_CASE(
    GlobalToBoundTrackParameters,
    surfaces* posSymmetric* posSymmetric* ts* phis* thetas* ps* qsNonZero,
    surface, l0, l1, time, phiInput, theta, p, q) {
  // phi is ill-defined in forward/backward tracks
  const auto phi = ((0 < theta) and (theta < M_PI)) ? phiInput : 0.0;
  const auto qOverP = q / p;

  GeometryContext geoCtx;
  Vector2D loc(l0, l1);
  Vector3D dir = makeDirectionUnitFromPhiTheta(phi, theta);
  // transform reference position
  Vector3D pos = surface->localToGlobal(geoCtx, loc, dir);

  // convert free parameters to bound parameters
  {
    BOOST_TEST_INFO("Transform free parameters vector onto surface "
                    << surface->name());

    FreeVector fv = FreeVector::Zero();
    fv[eFreePos0] = pos[ePos0];
    fv[eFreePos1] = pos[ePos1];
    fv[eFreePos2] = pos[ePos2];
    fv[eFreeTime] = time;
    fv[eFreeDir0] = dir[eMom0];
    fv[eFreeDir1] = dir[eMom1];
    fv[eFreeDir2] = dir[eMom2];
    fv[eFreeQOverP] = qOverP;
    BoundVector bv =
        detail::transformFreeToBoundParameters(fv, *surface, geoCtx);
    CHECK_CLOSE_OR_SMALL(bv[eBoundLoc0], l0, eps, eps);
    CHECK_CLOSE_OR_SMALL(bv[eBoundLoc1], l1, eps, eps);
    CHECK_CLOSE_OR_SMALL(bv[eBoundTime], time, eps, eps);
    CHECK_CLOSE_OR_SMALL(bv[eBoundPhi], phi, eps, eps);
    CHECK_CLOSE_OR_SMALL(bv[eBoundTheta], theta, eps, eps);
    CHECK_CLOSE_OR_SMALL(bv[eBoundQOverP], qOverP, eps, eps);
  }
  // convert separate components to bound parameters
  {
    BOOST_TEST_INFO("Transform free parameters components onto surface "
                    << surface->name());

    BoundVector bv = detail::transformFreeToBoundParameters(
        pos, time, dir, qOverP, *surface, geoCtx);
    CHECK_CLOSE_OR_SMALL(bv[eBoundLoc0], l0, eps, eps);
    CHECK_CLOSE_OR_SMALL(bv[eBoundLoc1], l1, eps, eps);
    CHECK_CLOSE_OR_SMALL(bv[eBoundTime], time, eps, eps);
    CHECK_CLOSE_OR_SMALL(bv[eBoundPhi], phi, eps, eps);
    CHECK_CLOSE_OR_SMALL(bv[eBoundTheta], theta, eps, eps);
    CHECK_CLOSE_OR_SMALL(bv[eBoundQOverP], qOverP, eps, eps);
  }
}

BOOST_DATA_TEST_CASE(GlobalToCurvilinearParameters,
                     ts* phis* thetas* ps* qsNonZero, time, phiInput, theta, p,
                     q) {
  // phi is ill-defined in forward/backward tracks
  const auto phi = ((0 < theta) and (theta < M_PI)) ? phiInput : 0.0;
  const auto qOverP = q / p;

  GeometryContext geoCtx;
  Vector3D dir = makeDirectionUnitFromPhiTheta(phi, theta);

  // convert w/ direction
  {
    BoundVector bv =
        detail::transformFreeToCurvilinearParameters(time, dir, qOverP);
    CHECK_SMALL(bv[eBoundLoc0], eps);
    CHECK_SMALL(bv[eBoundLoc1], eps);
    CHECK_CLOSE_OR_SMALL(bv[eBoundTime], time, eps, eps);
    CHECK_CLOSE_OR_SMALL(bv[eBoundPhi], phi, eps, eps);
    CHECK_CLOSE_OR_SMALL(bv[eBoundTheta], theta, eps, eps);
    CHECK_CLOSE_OR_SMALL(bv[eBoundQOverP], qOverP, eps, eps);
  }
  // convert w/ angles
  {
    BoundVector bv =
        detail::transformFreeToCurvilinearParameters(time, phi, theta, qOverP);
    CHECK_SMALL(bv[eBoundLoc0], eps);
    CHECK_SMALL(bv[eBoundLoc1], eps);
    CHECK_CLOSE_OR_SMALL(bv[eBoundTime], time, eps, eps);
    CHECK_CLOSE_OR_SMALL(bv[eBoundPhi], phi, eps, eps);
    CHECK_CLOSE_OR_SMALL(bv[eBoundTheta], theta, eps, eps);
    CHECK_CLOSE_OR_SMALL(bv[eBoundQOverP], qOverP, eps, eps);
  }
}

BOOST_AUTO_TEST_SUITE_END()

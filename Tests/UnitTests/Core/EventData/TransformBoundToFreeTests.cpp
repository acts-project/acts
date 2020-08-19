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

#include "Acts/EventData/detail/TransformationBoundToFree.hpp"
#include "Acts/Tests/CommonHelpers/FloatComparisons.hpp"
#include "Acts/Utilities/UnitVectors.hpp"
#include "Acts/Utilities/Units.hpp"
#include "TrackParametersTestData.hpp"

using namespace Acts;
using namespace Acts::UnitLiterals;

namespace {
constexpr auto eps = std::numeric_limits<FreeParametersScalar>::epsilon();
}

BOOST_AUTO_TEST_SUITE(TransformBoundToFree)

BOOST_DATA_TEST_CASE(
    Parameters, surfaces* posSymmetric* posSymmetric* ts* phis* thetas* qOverPs,
    surface, l0, l1, time, phi, theta, qOverP) {
  GeometryContext geoCtx;

  Vector2D loc(l0, l1);
  Vector3D pos = Vector3D::Zero();
  Vector3D dir = makeDirectionUnitFromPhiTheta(phi, theta);
  // transform reference position
  surface->localToGlobal(geoCtx, loc, dir, pos);

  // construct bound parameters
  BoundVector bv = BoundVector::Zero();
  bv[eBoundLoc0] = l0;
  bv[eBoundLoc1] = l1;
  bv[eBoundTime] = time;
  bv[eBoundPhi] = phi;
  bv[eBoundTheta] = theta;
  bv[eBoundQOverP] = qOverP;

  // convert to free parameters
  FreeVector fv = detail::transformBoundToFreeParameters(*surface, geoCtx, bv);

  CHECK_CLOSE_OR_SMALL(fv.segment<3>(eFreePos0), pos, eps, eps);
  CHECK_CLOSE_OR_SMALL(fv[eFreeTime], bv[eBoundTime], eps, eps);
  CHECK_CLOSE_REL(fv.segment<3>(eFreeDir0).norm(), 1, eps);
  CHECK_CLOSE_OR_SMALL(fv.segment<3>(eFreeDir0), dir, eps, eps);
  CHECK_CLOSE_OR_SMALL(fv[eFreeQOverP], bv[eBoundQOverP], eps, eps);
}

BOOST_AUTO_TEST_SUITE_END()

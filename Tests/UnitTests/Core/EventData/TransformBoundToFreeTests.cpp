// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/test/data/test_case.hpp>
#include <boost/test/unit_test.hpp>

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/Definitions/Units.hpp"
#include "Acts/EventData/detail/TransformationBoundToFree.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Tests/CommonHelpers/FloatComparisons.hpp"
#include "Acts/Utilities/UnitVectors.hpp"

#include <algorithm>
#include <limits>
#include <utility>
#include <vector>

#include "TrackParametersDatasets.hpp"

using namespace Acts;
using namespace Acts::UnitLiterals;

namespace {
constexpr auto eps = std::numeric_limits<ActsScalar>::epsilon();
}

BOOST_AUTO_TEST_SUITE(TransformBoundToFree)

BOOST_DATA_TEST_CASE(
    Parameters,
    surfaces* posSymmetric* posSymmetric* ts* phis* thetas* ps* qsNonZero,
    surface, l0, l1, time, phi, theta, p, q) {
  GeometryContext geoCtx;

  Vector2 loc(l0, l1);
  Vector3 dir = makeDirectionFromPhiTheta(phi, theta);
  // transform reference position
  Vector3 pos = surface->localToGlobal(geoCtx, loc, dir);

  const auto qOverP = q / p;

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

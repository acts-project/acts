// This file is part of the Acts project.
//
// Copyright (C) 2017-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/test/data/test_case.hpp>
#include <boost/test/unit_test.hpp>

#include "Acts/EventData/NeutralTrackParameters.hpp"
#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/Surfaces/ConeSurface.hpp"
#include "Acts/Surfaces/CylinderSurface.hpp"
#include "Acts/Surfaces/DiscSurface.hpp"
#include "Acts/Surfaces/PerigeeSurface.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/Surfaces/StrawSurface.hpp"
#include "Acts/Tests/CommonHelpers/FloatComparisons.hpp"
#include "Acts/Utilities/Units.hpp"
#include "Acts/Utilities/detail/periodic.hpp"

#include "TrackParametersTestData.hpp"

namespace {

namespace bdata = boost::unit_test::data;
using namespace Acts;
using namespace Acts::UnitLiterals;

static constexpr auto eps =
    8 * std::numeric_limits<BoundParametersScalar>::epsilon();
static const GeometryContext geoCtx;

template <typename charge_t>
void checkParameters(const SingleBoundTrackParameters<charge_t>& params,
                     double l0, double l1, double time, double phi,
                     double theta, double qOverP, const Vector3D& pos,
                     const Vector3D& mom, double q) {
  // native values
  CHECK_CLOSE_OR_SMALL(params.template get<eBoundLoc0>(), l0, eps, eps);
  CHECK_CLOSE_OR_SMALL(params.template get<eBoundLoc1>(), l1, eps, eps);
  CHECK_CLOSE_OR_SMALL(params.template get<eBoundTime>(), time, eps, eps);
  CHECK_CLOSE_OR_SMALL(detail::radian_sym(params.template get<eBoundPhi>()),
                       detail::radian_sym(phi), eps, eps);
  CHECK_CLOSE_OR_SMALL(params.template get<eBoundTheta>(), theta, eps, eps);
  CHECK_CLOSE_OR_SMALL(params.template get<eBoundQOverP>(), qOverP, eps, eps);
  // convenience accessors
  CHECK_CLOSE_OR_SMALL(params.position(geoCtx), pos, eps, eps);
  CHECK_CLOSE_OR_SMALL(params.time(), time, eps, eps);
  CHECK_CLOSE_OR_SMALL(params.momentum(), mom, eps, eps);
  BOOST_CHECK_EQUAL(params.charge(), q);
}

void runTest(std::shared_ptr<const Surface> surface, double l0, double l1,
             double time, double phi, double theta, double p) {
  // phi is ill-defined in forward/backward tracks
  phi = ((0 < theta) and (theta < M_PI)) ? phi : 0.0;

  // global direction/momentum for reference
  const Vector3D dir = makeDirectionUnitFromPhiTheta(phi, theta);
  const Vector3D mom = p * dir;

  // convert local-to-global for reference
  const Vector2D loc(l0, l1);
  Vector3D pos = Vector3D::Zero();
  surface->localToGlobal(geoCtx, loc, dir, pos);

  // positively charged from local vector
  {
    BoundVector vector = BoundVector::Zero();
    vector[eBoundLoc0] = l0;
    vector[eBoundLoc1] = l1;
    vector[eBoundTime] = time;
    vector[eBoundPhi] = phi;
    vector[eBoundTheta] = theta;
    vector[eBoundQOverP] = 1_e / p;
    BoundParameters params(surface, vector);
    checkParameters(params, l0, l1, time, phi, theta, 1_e / p, pos, mom, 1_e);
  }
  // positively charged from global information
  {
    BoundParameters params(geoCtx, std::nullopt, pos, mom, 1_e, time, surface);
    checkParameters(params, l0, l1, time, phi, theta, 1_e / p, pos, mom, 1_e);
  }
  // negatively charged from local vector
  {
    BoundVector vector = BoundVector::Zero();
    vector[eBoundLoc0] = l0;
    vector[eBoundLoc1] = l1;
    vector[eBoundTime] = time;
    vector[eBoundPhi] = phi;
    vector[eBoundTheta] = theta;
    vector[eBoundQOverP] = -1_e / p;
    BoundParameters params(surface, vector);
    checkParameters(params, l0, l1, time, phi, theta, -1_e / p, pos, mom, -1_e);
  }
  // negatively charged from global information
  {
    BoundParameters params(geoCtx, std::nullopt, pos, mom, -1_e, time, surface);
    checkParameters(params, l0, l1, time, phi, theta, -1_e / p, pos, mom, -1_e);
  }
  // neutral parameters from local vector
  {
    BoundVector vector = BoundVector::Zero();
    vector[eBoundLoc0] = l0;
    vector[eBoundLoc1] = l1;
    vector[eBoundTime] = time;
    vector[eBoundPhi] = phi;
    vector[eBoundTheta] = theta;
    vector[eBoundQOverP] = 1 / p;
    NeutralBoundTrackParameters params(surface, vector);
    checkParameters(params, l0, l1, time, phi, theta, 1 / p, pos, mom, 0);
  }
  // neutral parameters from global information
  {
    NeutralBoundTrackParameters params(geoCtx, std::nullopt, pos, mom, time,
                                       surface);
    checkParameters(params, l0, l1, time, phi, theta, 1 / p, pos, mom, 0);
  }
}

std::shared_ptr<Transform3D> makeTransformIdentity() {
  return std::make_shared<Transform3D>(Transform3D::Identity());
}

// different surfaces
// parameters must be chosen such that all possible local positions (as defined
// in the datasets header) represent valid points on the surface.
const auto cones = bdata::make({
    Surface::makeShared<ConeSurface>(makeTransformIdentity(),
                                     0.5 /* opening angle */),
});
const auto cylinders = bdata::make({
    Surface::makeShared<CylinderSurface>(makeTransformIdentity(),
                                         10.0 /* radius */, 100 /* half z */),
});
const auto discs = bdata::make({
    Surface::makeShared<DiscSurface>(makeTransformIdentity(),
                                     0 /* radius min */, 100 /* radius max */),
});
const auto perigees = bdata::make({
    Surface::makeShared<PerigeeSurface>(Vector3D(0, 0, -1.5)),
});
const auto planes = bdata::make({
    Surface::makeShared<PlaneSurface>(Vector3D(1, 2, 3), Vector3D::UnitX()),
    Surface::makeShared<PlaneSurface>(Vector3D(-2, -3, -4), Vector3D::UnitY()),
    Surface::makeShared<PlaneSurface>(Vector3D(3, -4, 5), Vector3D::UnitZ()),
});
const auto straws = bdata::make({
    Surface::makeShared<StrawSurface>(makeTransformIdentity(), 2.0 /* radius */,
                                      200.0 /* half z */),
});

}  // namespace

BOOST_AUTO_TEST_SUITE(BoundTrackParameters)

BOOST_DATA_TEST_CASE(ConeSurface,
                     cones* posAngle* posPositiveNonzero* ts* phis* thetas* ps,
                     surface, lphi, lz, time, phi, theta, p) {
  // TODO extend lz to zero after fixing the transform implementation
  // local parameter r*phi has limits that depend on the z position
  const auto r = lz * surface->bounds().tanAlpha();
  // local coordinates are singular at z = 0 -> normalize local r*phi
  runTest(surface, (0 < lz) ? (r * lphi) : 0.0, lz, time, phi, theta, p);
}

BOOST_DATA_TEST_CASE(
    CylinderSurface,
    cylinders* posSymmetric* posSymmetric* ts* phis* thetas* ps, surface, lrphi,
    lz, time, phi, theta, p) {
  runTest(surface, lrphi, lz, time, phi, theta, p);
}

BOOST_DATA_TEST_CASE(DiscSurface,
                     discs* posPositive* posAngle* ts* phis* thetas* ps,
                     surface, lr, lphi, time, phi, theta, p) {
  // local coordinates are singular at r = 0 -> normalize local phi
  runTest(surface, lr, (0 < lr) ? lphi : 0.0, time, phi, theta, p);
}

BOOST_DATA_TEST_CASE(
    PerigeeSurface,
    perigees* posSymmetric* posSymmetric* ts* phis* thetasNoForwardBackward* ps,
    surface, d0, z0, time, phi, theta, p) {
  // TODO extend theta to forward/back extreme cases fixing the transform
  runTest(surface, d0, z0, time, phi, theta, p);
}

BOOST_DATA_TEST_CASE(PlaneSurface,
                     planes* posSymmetric* posSymmetric* ts* phis* thetas* ps,
                     surface, l0, l1, time, phi, theta, p) {
  runTest(surface, l0, l1, time, phi, theta, p);
}

BOOST_DATA_TEST_CASE(
    StrawSurface,
    straws* posPositive* posSymmetric* ts* phis* thetasNoForwardBackward* ps,
    surface, lr, lz, time, phi, theta, p) {
  // TODO extend theta to forward/back extreme cases fixing the transform
  runTest(surface, lr, lz, time, phi, theta, p);
}

BOOST_AUTO_TEST_SUITE_END()

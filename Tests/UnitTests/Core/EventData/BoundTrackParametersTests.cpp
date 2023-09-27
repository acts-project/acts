// This file is part of the Acts project.
//
// Copyright (C) 2017-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/test/data/test_case.hpp>
#include <boost/test/unit_test.hpp>

#include "Acts/Definitions/Units.hpp"
#include "Acts/EventData/NeutralTrackParameters.hpp"
#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/Surfaces/ConeSurface.hpp"
#include "Acts/Surfaces/CylinderSurface.hpp"
#include "Acts/Surfaces/DiscSurface.hpp"
#include "Acts/Surfaces/PerigeeSurface.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/Surfaces/StrawSurface.hpp"
#include "Acts/Tests/CommonHelpers/FloatComparisons.hpp"
#include "Acts/Utilities/Helpers.hpp"
#include "Acts/Utilities/UnitVectors.hpp"
#include "Acts/Utilities/detail/periodic.hpp"

#include "TrackParametersDatasets.hpp"

namespace {

namespace bdata = boost::unit_test::data;
using namespace Acts;
using namespace Acts::UnitLiterals;
using AnyBoundTrackParameters = SingleBoundTrackParameters<AnyCharge>;

constexpr auto eps = 8 * std::numeric_limits<ActsScalar>::epsilon();
const GeometryContext geoCtx;
const BoundSymMatrix cov = BoundSymMatrix::Identity();

template <typename charge_t>
void checkParameters(const SingleBoundTrackParameters<charge_t>& params,
                     double l0, double l1, double time, double phi,
                     double theta, double p, double q, const Vector3& pos,
                     const Vector3& unitDir) {
  const auto qOverP = (q != 0) ? (q / p) : (1 / p);
  const auto pos4 = VectorHelpers::makeVector4(pos, time);

  // native values
  CHECK_CLOSE_OR_SMALL(params.template get<eBoundLoc0>(), l0, eps, eps);
  CHECK_CLOSE_OR_SMALL(params.template get<eBoundLoc1>(), l1, eps, eps);
  CHECK_CLOSE_OR_SMALL(params.template get<eBoundTime>(), time, eps, eps);
  CHECK_CLOSE_OR_SMALL(detail::radian_sym(params.template get<eBoundPhi>()),
                       detail::radian_sym(phi), eps, eps);
  CHECK_CLOSE_OR_SMALL(params.template get<eBoundTheta>(), theta, eps, eps);
  CHECK_CLOSE_OR_SMALL(params.template get<eBoundQOverP>(), qOverP, eps, eps);
  // convenience accessors
  CHECK_CLOSE_OR_SMALL(params.fourPosition(geoCtx), pos4, eps, eps);
  CHECK_CLOSE_OR_SMALL(params.position(geoCtx), pos, eps, eps);
  CHECK_CLOSE_OR_SMALL(params.time(), time, eps, eps);
  CHECK_CLOSE_OR_SMALL(params.unitDirection(), unitDir, eps, eps);
  CHECK_CLOSE_OR_SMALL(params.absoluteMomentum(), p, eps, eps);
  CHECK_CLOSE_OR_SMALL(params.transverseMomentum(), p * std::sin(theta), eps,
                       eps);
  CHECK_CLOSE_OR_SMALL(params.momentum(), p * unitDir, eps, eps);
  BOOST_CHECK_EQUAL(params.charge(), q);
}

void runTest(const std::shared_ptr<const Surface>& surface, double l0,
             double l1, double time, double phi, double theta, double p) {
  // phi is ill-defined in forward/backward tracks
  phi = ((0 < theta) and (theta < M_PI)) ? phi : 0.0;

  // global direction for reference
  const Vector3 dir = makeDirectionUnitFromPhiTheta(phi, theta);
  // convert local-to-global for reference
  const Vector2 loc(l0, l1);
  const Vector3 pos = surface->localToGlobal(geoCtx, loc, dir);
  // global four-position as input
  Vector4 pos4;
  pos4.segment<3>(ePos0) = pos;
  pos4[eTime] = time;

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
    checkParameters(params, l0, l1, time, phi, theta, p, 0_e, pos, dir);
    BOOST_CHECK(not params.covariance());

    // reassign w/ covariance
    params = NeutralBoundTrackParameters(surface, vector, cov);
    checkParameters(params, l0, l1, time, phi, theta, p, 0_e, pos, dir);
    BOOST_CHECK(params.covariance());
    BOOST_CHECK_EQUAL(params.covariance().value(), cov);
  }
  // negative charged parameters from local vector
  {
    BoundVector vector = BoundVector::Zero();
    vector[eBoundLoc0] = l0;
    vector[eBoundLoc1] = l1;
    vector[eBoundTime] = time;
    vector[eBoundPhi] = phi;
    vector[eBoundTheta] = theta;
    vector[eBoundQOverP] = -1_e / p;
    BoundTrackParameters params(surface, vector);
    checkParameters(params, l0, l1, time, phi, theta, p, -1_e, pos, dir);
    BOOST_CHECK(not params.covariance());

    // reassign w/ covariance
    params = BoundTrackParameters(surface, vector, cov);
    checkParameters(params, l0, l1, time, phi, theta, p, -1_e, pos, dir);
    BOOST_CHECK(params.covariance());
    BOOST_CHECK_EQUAL(params.covariance().value(), cov);
  }
  // positive charged parameters from local vector
  {
    BoundVector vector = BoundVector::Zero();
    vector[eBoundLoc0] = l0;
    vector[eBoundLoc1] = l1;
    vector[eBoundTime] = time;
    vector[eBoundPhi] = phi;
    vector[eBoundTheta] = theta;
    vector[eBoundQOverP] = 1_e / p;
    BoundTrackParameters params(surface, vector);
    checkParameters(params, l0, l1, time, phi, theta, p, 1_e, pos, dir);
    BOOST_CHECK(not params.covariance());

    // reassign w/ covariance
    params = BoundTrackParameters(surface, vector, cov);
    checkParameters(params, l0, l1, time, phi, theta, p, 1_e, pos, dir);
    BOOST_CHECK(params.covariance());
    BOOST_CHECK_EQUAL(params.covariance().value(), cov);
  }
  // double-negative charged any parameters from local vector
  {
    BoundVector vector = BoundVector::Zero();
    vector[eBoundLoc0] = l0;
    vector[eBoundLoc1] = l1;
    vector[eBoundTime] = time;
    vector[eBoundPhi] = phi;
    vector[eBoundTheta] = theta;
    vector[eBoundQOverP] = -2_e / p;
    AnyBoundTrackParameters params(surface, vector, -2_e);
    checkParameters(params, l0, l1, time, phi, theta, p, -2_e, pos, dir);
    BOOST_CHECK(not params.covariance());

    // reassign w/ covariance
    params = AnyBoundTrackParameters(surface, vector, -2_e, cov);
    checkParameters(params, l0, l1, time, phi, theta, p, -2_e, pos, dir);
    BOOST_CHECK(params.covariance());
    BOOST_CHECK_EQUAL(params.covariance().value(), cov);
  }
  // neutral parameters from global information
  {
    auto params =
        NeutralBoundTrackParameters::create(surface, geoCtx, pos4, dir, 1 / p)
            .value();
    checkParameters(params, l0, l1, time, phi, theta, p, 0_e, pos, dir);
    BOOST_CHECK(not params.covariance());
  }
  // negative charged parameters from global information
  {
    auto params =
        BoundTrackParameters::create(surface, geoCtx, pos4, dir, -1_e / p)
            .value();
    checkParameters(params, l0, l1, time, phi, theta, p, -1_e, pos, dir);
    BOOST_CHECK(not params.covariance());
  }
  // positive charged parameters from global information
  {
    auto params =
        BoundTrackParameters::create(surface, geoCtx, pos4, dir, 1_e / p)
            .value();
    checkParameters(params, l0, l1, time, phi, theta, p, 1_e, pos, dir);
    BOOST_CHECK(not params.covariance());
  }
  // neutral any parameters from global information
  {
    auto params =
        AnyBoundTrackParameters::create(surface, geoCtx, pos4, dir, p, 0_e)
            .value();
    checkParameters(params, l0, l1, time, phi, theta, p, 0_e, pos, dir);
    BOOST_CHECK(not params.covariance());
  }
  // double-negative any parameters from global information
  {
    auto params =
        AnyBoundTrackParameters::create(surface, geoCtx, pos4, dir, p, -2_e)
            .value();
    checkParameters(params, l0, l1, time, phi, theta, p, -2_e, pos, dir);
    BOOST_CHECK(not params.covariance());
  }
  // triple-positive any parameters from global information
  {
    auto params =
        AnyBoundTrackParameters::create(surface, geoCtx, pos4, dir, p, 3_e)
            .value();
    checkParameters(params, l0, l1, time, phi, theta, p, 3_e, pos, dir);
    BOOST_CHECK(not params.covariance());
  }
}

// different surfaces
// parameters must be chosen such that all possible local positions (as defined
// in the datasets header) represent valid points on the surface.
const auto cones = bdata::make({
    Surface::makeShared<ConeSurface>(Transform3::Identity(),
                                     0.5 /* opening angle */),
});
const auto cylinders = bdata::make({
    Surface::makeShared<CylinderSurface>(Transform3::Identity(),
                                         10.0 /* radius */, 100 /* half z */),
});
const auto discs = bdata::make({
    Surface::makeShared<DiscSurface>(Transform3::Identity(), 0 /* radius min */,
                                     100 /* radius max */),
});
const auto perigees = bdata::make({
    Surface::makeShared<PerigeeSurface>(Vector3(0, 0, -1.5)),
});
const auto planes = bdata::make({
    Surface::makeShared<PlaneSurface>(Vector3(1, 2, 3), Vector3::UnitX()),
    Surface::makeShared<PlaneSurface>(Vector3(-2, -3, -4), Vector3::UnitY()),
    Surface::makeShared<PlaneSurface>(Vector3(3, -4, 5), Vector3::UnitZ()),
});
const auto straws = bdata::make({
    Surface::makeShared<StrawSurface>(Transform3::Identity(), 2.0 /* radius */,
                                      200.0 /* half z */),
});

}  // namespace

BOOST_AUTO_TEST_SUITE(EventDataBoundTrackParameters)

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

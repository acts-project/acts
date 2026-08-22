// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/data/test_case.hpp>
#include <boost/test/unit_test.hpp>

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/Definitions/Units.hpp"
#include "Acts/EventData/BoundTrackParameters.hpp"
#include "Acts/EventData/ParticleHypothesis.hpp"
#include "Acts/EventData/TransformationHelpers.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Utilities/UnitVectors.hpp"
#include "ActsTests/CommonHelpers/FloatComparisons.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <numbers>
#include <optional>
#include <utility>
#include <vector>

#include "TrackParametersDatasets.hpp"

using namespace Acts;
using namespace Acts::UnitLiterals;

namespace {
constexpr auto eps = std::numeric_limits<double>::epsilon();
}

namespace ActsTests {

BOOST_AUTO_TEST_SUITE(EventDataSuite)

BOOST_DATA_TEST_CASE(Parameters,
                     surfaces * posSymmetric * posSymmetric * ts * phis *
                         thetas * ps * qsNonZero,
                     surface, l0, l1, time, phi, theta, p, q) {
  auto geoCtx = GeometryContext::dangerouslyDefaultConstruct();

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
  FreeVector fv = transformBoundToFreeParameters(*surface, geoCtx, bv);

  CHECK_CLOSE_OR_SMALL(fv.segment<3>(eFreePos0), pos, eps, eps);
  CHECK_CLOSE_OR_SMALL(fv[eFreeTime], bv[eBoundTime], eps, eps);
  CHECK_CLOSE_REL(fv.segment<3>(eFreeDir0).norm(), 1, eps);
  CHECK_CLOSE_OR_SMALL(fv.segment<3>(eFreeDir0), dir, eps, eps);
  CHECK_CLOSE_OR_SMALL(fv[eFreeQOverP], bv[eBoundQOverP], eps, eps);
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(TransformFreeToBound)

BOOST_DATA_TEST_CASE(GlobalToBoundTrackParameters,
                     surfaces * posSymmetric * posSymmetric * ts * phis *
                         thetas * ps * qsNonZero,
                     surface, l0, l1, time, phiInput, theta, p, q) {
  // phi is ill-defined in forward/backward tracks
  const auto phi = ((0 < theta) && (theta < std::numbers::pi)) ? phiInput : 0.;
  const auto qOverP = q / p;

  auto geoCtx = GeometryContext::dangerouslyDefaultConstruct();
  Vector2 loc(l0, l1);
  Vector3 dir = makeDirectionFromPhiTheta(phi, theta);
  // transform reference position
  Vector3 pos = surface->localToGlobal(geoCtx, loc, dir);

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
        transformFreeToBoundParameters(fv, *surface, geoCtx).value();
    CHECK_CLOSE_OR_SMALL(bv[eBoundLoc0], l0, eps, eps);
    CHECK_CLOSE_OR_SMALL(bv[eBoundLoc1], l1, eps, eps);
    CHECK_CLOSE_OR_SMALL(bv[eBoundTime], time, eps, eps);
    CHECK_CLOSE_OR_SMALL(bv[eBoundPhi], phi, eps, eps);
    CHECK_CLOSE_OR_SMALL(bv[eBoundTheta], theta, eps, eps);
    CHECK_CLOSE_OR_SMALL(bv[eBoundQOverP], qOverP, eps, eps);
  }

  // Assert failure when trying to convert a position that is not on-surface.
  {
    Vector3 posOff = pos + surface->normal(geoCtx, loc) * 0.5;
    BOOST_TEST_INFO("Transform free parameters vector onto surface "
                    << surface->name());

    FreeVector fv = FreeVector::Zero();
    fv[eFreePos0] = posOff[ePos0];
    fv[eFreePos1] = posOff[ePos1];
    fv[eFreePos2] = posOff[ePos2];
    fv[eFreeTime] = time;
    fv[eFreeDir0] = dir[eMom0];
    fv[eFreeDir1] = dir[eMom1];
    fv[eFreeDir2] = dir[eMom2];
    fv[eFreeQOverP] = qOverP;
    auto res = transformFreeToBoundParameters(fv, *surface, geoCtx);
    BOOST_CHECK(!res.ok());
  }

  // convert separate components to bound parameters
  {
    BOOST_TEST_INFO("Transform free parameters components onto surface "
                    << surface->name());

    BoundVector bv =
        transformFreeToBoundParameters(pos, time, dir, qOverP, *surface, geoCtx)
            .value();
    CHECK_CLOSE_OR_SMALL(bv[eBoundLoc0], l0, eps, eps);
    CHECK_CLOSE_OR_SMALL(bv[eBoundLoc1], l1, eps, eps);
    CHECK_CLOSE_OR_SMALL(bv[eBoundTime], time, eps, eps);
    CHECK_CLOSE_OR_SMALL(bv[eBoundPhi], phi, eps, eps);
    CHECK_CLOSE_OR_SMALL(bv[eBoundTheta], theta, eps, eps);
    CHECK_CLOSE_OR_SMALL(bv[eBoundQOverP], qOverP, eps, eps);
  }

  // Assert failure when trying to convert a position that is not on-surface.
  {
    BOOST_TEST_INFO("Transform free parameters components onto surface "
                    << surface->name());

    Vector3 posOff = pos + surface->normal(geoCtx, loc) * 0.5;
    auto res = transformFreeToBoundParameters(posOff, time, dir, qOverP,
                                              *surface, geoCtx);
    BOOST_CHECK(!res.ok());
  }
}

BOOST_DATA_TEST_CASE(GlobalToCurvilinearParameters,
                     ts * phis * thetas * ps * qsNonZero, time, phiInput, theta,
                     p, q) {
  // phi is ill-defined in forward/backward tracks
  const auto phi = ((0 < theta) && (theta < std::numbers::pi)) ? phiInput : 0.;
  const auto qOverP = q / p;

  auto geoCtx = GeometryContext::dangerouslyDefaultConstruct();
  Vector3 dir = makeDirectionFromPhiTheta(phi, theta);

  // convert w/ direction
  {
    BoundVector bv = transformFreeToCurvilinearParameters(time, dir, qOverP);
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
        transformFreeToCurvilinearParameters(time, phi, theta, qOverP);
    CHECK_SMALL(bv[eBoundLoc0], eps);
    CHECK_SMALL(bv[eBoundLoc1], eps);
    CHECK_CLOSE_OR_SMALL(bv[eBoundTime], time, eps, eps);
    CHECK_CLOSE_OR_SMALL(bv[eBoundPhi], phi, eps, eps);
    CHECK_CLOSE_OR_SMALL(bv[eBoundTheta], theta, eps, eps);
    CHECK_CLOSE_OR_SMALL(bv[eBoundQOverP], qOverP, eps, eps);
  }
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(TransformBoundToCartesian)

BOOST_DATA_TEST_CASE(BoundToCartesianFourPositionMomentum,
                     surfaces * posSymmetric * posSymmetric * ts * phis *
                         thetasNoForwardBackward * ps * qsNonZero,
                     surface, l0, l1, time, phi, theta, p, q) {
  // Forward/backward tracks (theta = 0, pi) are excluded: phi is ill-defined
  // there and the bound-to-free direction Jacobian is singular (division by
  // sin(theta) in sphericalToFreeDirectionJacobian), so
  // transformBoundToCartesianFourPositionMomentum cannot be evaluated.
  const auto qOverP = q / p;

  auto geoCtx = GeometryContext::dangerouslyDefaultConstruct();
  Vector2 loc(l0, l1);
  Vector3 dir = makeDirectionFromPhiTheta(phi, theta);
  Vector3 pos = surface->localToGlobal(geoCtx, loc, dir);

  BoundVector bv = BoundVector::Zero();
  bv[eBoundLoc0] = l0;
  bv[eBoundLoc1] = l1;
  bv[eBoundTime] = time;
  bv[eBoundPhi] = phi;
  bv[eBoundTheta] = theta;
  bv[eBoundQOverP] = qOverP;

  BoundTrackParameters params(surface, bv, BoundMatrix::Identity(),
                              ParticleHypothesis::pionLike(std::abs(q)));

  Vector3 momentum;
  auto [pos4, cov7] =
      transformBoundToCartesianFourPositionMomentum(geoCtx, params, momentum);

  CHECK_CLOSE_OR_SMALL(pos4.segment<3>(ePos0), pos, eps, eps);
  CHECK_CLOSE_OR_SMALL(pos4[eTime], time, eps, eps);
  // "small" scales with p: rounding noise in near-zero components of dir
  // (e.g. sin(pi) != 0 exactly) gets amplified by the momentum magnitude.
  CHECK_CLOSE_OR_SMALL(momentum, p * dir, eps, p * eps);

  // the propagated covariance must remain symmetric; the absolute floor
  // scales with the matrix magnitude for the same reason as above.
  CHECK_CLOSE_OR_SMALL(cov7, cov7.transpose(), eps,
                       cov7.cwiseAbs().maxCoeff() * eps);
}

BOOST_DATA_TEST_CASE(BoundToCartesianCovarianceConsistency,
                     surfaces * phis * thetasNoForwardBackward * ps * qsNonZero,
                     surface, phi, theta, p, q) {
  auto geoCtx = GeometryContext::dangerouslyDefaultConstruct();
  const auto qOverP = q / p;
  const auto hypothesis = ParticleHypothesis::pionLike(std::abs(q));

  BoundVector bv = BoundVector::Zero();
  bv[eBoundLoc0] = 0.3;
  bv[eBoundLoc1] = -0.2;
  bv[eBoundTime] = 1.0;
  bv[eBoundPhi] = phi;
  bv[eBoundTheta] = theta;
  bv[eBoundQOverP] = qOverP;

  const BoundMatrix boundCov = BoundMatrix::Identity();
  BoundTrackParameters params(surface, bv, boundCov, hypothesis);

  Vector3 momentum;
  auto [pos4, cov7] =
      transformBoundToCartesianFourPositionMomentum(geoCtx, params, momentum);

  // Build a numerical Jacobian d(x,y,z,t,px,py,pz)/d(BoundVector) with
  // central differences and cross-check the analytic covariance transport.
  Matrix<7, 6> numJac = Matrix<7, 6>::Zero();
  for (int i = 0; i < 6; ++i) {
    const double h = 1e-6 * std::max(1.0, std::abs(bv[i]));

    BoundVector bvPlus = bv;
    bvPlus[i] += h;
    BoundVector bvMinus = bv;
    bvMinus[i] -= h;

    BoundTrackParameters paramsPlus(surface, bvPlus, std::nullopt, hypothesis);
    BoundTrackParameters paramsMinus(surface, bvMinus, std::nullopt,
                                     hypothesis);

    numJac.block<4, 1>(0, i) =
        (paramsPlus.fourPosition(geoCtx) - paramsMinus.fourPosition(geoCtx)) /
        (2 * h);
    numJac.block<3, 1>(4, i) =
        (paramsPlus.momentum() - paramsMinus.momentum()) / (2 * h);
  }

  const SquareMatrix<7> cov7Numerical = numJac * boundCov * numJac.transpose();

  // Element-wise comparison akin to CHECK_CLOSE_COVARIANCE, but with an
  // absolute floor: some surface/direction combinations have a Cartesian
  // coordinate with exactly zero variance (e.g. along a plane surface's fixed
  // normal), which makes a purely relative tolerance ill-defined.
  for (int col = 0; col < 7; ++col) {
    for (int row = col; row < 7; ++row) {
      const double orderOfMagnitude =
          std::sqrt(cov7Numerical(row, row) * cov7Numerical(col, col));
      const double diff = std::abs(cov7(row, col) - cov7Numerical(row, col));
      BOOST_CHECK_SMALL(diff, std::max(1e-6 * orderOfMagnitude, 1e-9));
    }
  }
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace ActsTests

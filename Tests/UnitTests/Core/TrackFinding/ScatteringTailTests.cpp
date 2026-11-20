// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/TrackFinding/ScatteringTail.hpp"
#include "Acts/TrackFinding/detail/ScatteringTailState.hpp"
#include "ActsTests/CommonHelpers/FloatComparisons.hpp"

#include <cmath>
#include <numbers>

using namespace Acts;

namespace ActsTests {

BOOST_AUTO_TEST_SUITE(TrackFindingSuite)

BOOST_AUTO_TEST_CASE(ChiSquareIsMixtureLikelihood) {
  const Vector2 residual(0.3, -0.1);
  SquareMatrix2 core;
  core << 0.01, 0.002, 0.002, 0.02;
  SquareMatrix2 tail;
  tail << 0.2, 0.05, 0.05, 0.1;
  const double weight = 0.02;

  const auto gauss = [&](const SquareMatrix2& cov) {
    return std::exp(-0.5 * residual.dot(cov.inverse() * residual)) /
           (2 * std::numbers::pi * std::sqrt(cov.determinant()));
  };
  const double coreLikelihood = (1 - weight) * gauss(core);
  const double tailLikelihood = weight * gauss(core + tail);
  const double chi2Core = residual.dot(core.inverse() * residual);

  const auto [chi2, a] = scatteringTailChi2<2>(residual, core, tail,
                                               std::log(weight / (1 - weight)));
  CHECK_CLOSE_REL(chi2,
                  chi2Core - 2 * std::log((coreLikelihood + tailLikelihood) /
                                          coreLikelihood),
                  1e-10);
  CHECK_CLOSE_REL(1 / (1 + std::exp(-a)),
                  tailLikelihood / (coreLikelihood + tailLikelihood), 1e-10);
  BOOST_CHECK_LT(chi2, chi2Core);
}

BOOST_AUTO_TEST_CASE(TailStateAccumulatesAndTransports) {
  detail::ScatteringTailState state;
  const BoundMatrix predicted = BoundMatrix::Identity();
  BOOST_CHECK(!state.prior(predicted, 0).has_value());

  BoundMatrix jacobian = BoundMatrix::Identity();
  jacobian(eBoundLoc0, eBoundPhi) = 10;
  jacobian(eBoundLoc1, eBoundTheta) = 20;

  const ScatteringTailComponent tail{.weight = 0.02, .widthScale = 5};
  BoundMatrix expected = BoundMatrix::Zero();
  double noTail = 1;
  for (int i = 1; i <= 3; ++i) {
    state.transport(jacobian);
    expected = jacobian * expected * jacobian.transpose();
    state.add(tail, 1e-6 * i, 2e-6 * i);
    expected(eBoundPhi, eBoundPhi) += 0.02 * 24 * 1e-6 * i;
    expected(eBoundTheta, eBoundTheta) += 0.02 * 24 * 2e-6 * i;
    noTail *= 1 - tail.weight;
  }
  state.transport(jacobian);
  expected = jacobian * expected * jacobian.transpose() / (3 * 0.02);

  const auto prior = state.prior(predicted, 0);
  BOOST_REQUIRE(prior.has_value());
  CHECK_CLOSE_REL(prior->weight, 1 - noTail, 1e-12);
  BOOST_CHECK(prior->covariance.isApprox(expected, 1e-12));

  // the gate drops a tail that is small against the prediction
  BOOST_CHECK(!state.prior(1e3 * predicted, 0.05).has_value());

  state.reset();
  BOOST_CHECK(!state.prior(predicted, 0).has_value());
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace ActsTests

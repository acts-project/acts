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
#include "Acts/EventData/BoundTrackParameters.hpp"
#include "Acts/EventData/MultiComponentTrackParameters.hpp"
#include "Acts/EventData/ParticleHypothesis.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Surfaces/CurvilinearSurface.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/TrackFitting/GsfOptions.hpp"

#include <cmath>
#include <initializer_list>
#include <memory>
#include <numbers>
#include <optional>
#include <tuple>
#include <vector>

using namespace Acts;
const GeometryContext gctx = GeometryContext::dangerouslyDefaultConstruct();

static const auto particleHypothesis = ParticleHypothesis::pion();

namespace ActsTests {

BOOST_AUTO_TEST_SUITE(EventDataSuite)

BOOST_AUTO_TEST_CASE(test_constructors) {
  std::vector<std::tuple<double, BoundVector, BoundMatrix>> a;
  a.emplace_back(1.0, BoundVector::Ones(), BoundMatrix::Identity());

  std::vector<std::tuple<double, BoundVector, std::optional<BoundMatrix>>> b;
  b.emplace_back(1.0, BoundVector::Ones(), BoundMatrix::Identity());

  std::shared_ptr<PlaneSurface> surface =
      CurvilinearSurface(Vector3::Ones(), Vector3::Ones().normalized())
          .planeSurface();

  const auto ap = MultiComponentBoundTrackParameters(
      surface, a,
      [](const auto &cmp) {
        return std::tie(std::get<0>(cmp), std::get<1>(cmp), std::get<2>(cmp));
      },
      particleHypothesis);
  const auto bp = MultiComponentBoundTrackParameters(
      surface, b,
      [](const auto &cmp) {
        return std::tie(std::get<0>(cmp), std::get<1>(cmp), *std::get<2>(cmp));
      },
      particleHypothesis);
  const auto aps = MultiComponentBoundTrackParameters(
      surface, std::get<1>(a.front()), std::get<2>(a.front()),
      particleHypothesis);
  const auto bps = MultiComponentBoundTrackParameters(
      surface, std::get<1>(b.front()), std::get<2>(b.front()),
      particleHypothesis);

  BOOST_CHECK(b == ap.toComponents());
  BOOST_CHECK(ap.toComponents() == bp.toComponents());
  BOOST_CHECK(bp.toComponents() == aps.toComponents());
  BOOST_CHECK(aps.toComponents() == bps.toComponents());

  BOOST_CHECK_THROW(
      MultiComponentBoundTrackParameters(
          surface, std::vector<std::tuple<double, BoundVector, BoundMatrix>>{},
          [](const auto &cmp) {
            return std::tie(std::get<0>(cmp), std::get<1>(cmp),
                            std::get<2>(cmp));
          },
          particleHypothesis),
      std::invalid_argument);
}

BOOST_AUTO_TEST_CASE(test_accessors) {
  using cov_t = std::optional<BoundMatrix>;
  for (const auto &cov : {cov_t{}, cov_t{BoundMatrix::Identity()},
                          cov_t{BoundMatrix::Identity()}}) {
    std::shared_ptr<PlaneSurface> surface =
        CurvilinearSurface(Vector3::Ones(), Vector3::Ones().normalized())
            .planeSurface();

    const BoundTrackParameters single_pars(surface, BoundVector::Ones(), cov,
                                           particleHypothesis);

    const auto multi_pars = [&]() {
      MultiComponentBoundTrackParameters result(
          surface, single_pars.covariance().has_value(),
          single_pars.particleHypothesis());
      for (int i = 0; i < 4; ++i) {
        result.pushComponent(0.25, single_pars.parameters(),
                             single_pars.covariance());
      }
      return result;
    }();
    const auto multi_pars_as_single =
        multi_pars.merge(gctx, ComponentMergeMethod::eMean);

    BOOST_CHECK_EQUAL(multi_pars_as_single.absoluteMomentum(),
                      single_pars.absoluteMomentum());
    BOOST_CHECK_EQUAL(multi_pars_as_single.charge(), single_pars.charge());
    BOOST_CHECK_EQUAL(multi_pars_as_single.fourPosition(
                          GeometryContext::dangerouslyDefaultConstruct()),
                      single_pars.fourPosition(
                          GeometryContext::dangerouslyDefaultConstruct()));
    BOOST_CHECK_EQUAL(multi_pars_as_single.momentum(), single_pars.momentum());
    BOOST_CHECK_EQUAL(multi_pars_as_single.parameters(),
                      single_pars.parameters());
    BOOST_CHECK_EQUAL(
        multi_pars_as_single.position(
            GeometryContext::dangerouslyDefaultConstruct()),
        single_pars.position(GeometryContext::dangerouslyDefaultConstruct()));
    BOOST_CHECK_EQUAL(multi_pars_as_single.transverseMomentum(),
                      single_pars.transverseMomentum());
    BOOST_CHECK_EQUAL(multi_pars_as_single.direction(),
                      single_pars.direction());

    // Check the behaviour for std::nullopt or zero covariance
    if (cov.has_value()) {
      BOOST_CHECK_EQUAL(*multi_pars_as_single.covariance(),
                        *single_pars.covariance());
    } else {
      BOOST_CHECK(!multi_pars_as_single.covariance());
    }
  }
}
BOOST_AUTO_TEST_CASE(MergeAcrossPhiSeamWithAndWithoutCovariance) {
  const auto surface =
      Surface::makeShared<PlaneSurface>(Transform3::Identity());
  BoundVector left = BoundVector::Zero();
  left[eBoundLoc0] = -.3;
  left[eBoundDir0] = std::numbers::pi - .1;
  left[eBoundDir1] = 1.;
  left[eBoundQOverP] = .1;
  BoundVector right = left;
  right[eBoundLoc0] = .3;
  right[eBoundDir0] = -std::numbers::pi + .1;
  const BoundMatrix componentCovariance = .04 * BoundMatrix::Identity();

  for (bool withCovariance : {false, true}) {
    MultiComponentBoundTrackParameters mixture(surface, withCovariance,
                                               particleHypothesis);
    const std::optional<BoundMatrix> covariance =
        withCovariance ? std::optional{componentCovariance} : std::nullopt;
    mixture.pushComponent(.5, left, covariance);
    mixture.pushComponent(.5, right, covariance);
    const auto merged = mixture.merge(gctx, ComponentMergeMethod::eMean);
    BOOST_CHECK_SMALL(
        std::abs(merged.parameters()[eBoundDir0]) - std::numbers::pi, 1e-12);
    BOOST_CHECK_SMALL(merged.parameters()[eBoundLoc0], 1e-12);
    BOOST_CHECK_EQUAL(merged.covariance().has_value(), withCovariance);
    if (withCovariance) {
      BOOST_CHECK_SMALL((*merged.covariance())(eBoundDir0, eBoundDir0) - .05,
                        1e-12);
      BOOST_CHECK_SMALL((*merged.covariance())(eBoundLoc0, eBoundDir0) - .03,
                        1e-12);
    }

    mixture.clear();
    mixture.pushComponent(.25, left, covariance);
    mixture.pushComponent(.75, right, covariance);
    const auto mode = mixture.merge(gctx, ComponentMergeMethod::eMaxWeight);
    BOOST_CHECK(mode.parameters().isApprox(right, 1e-12));
  }
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace ActsTests

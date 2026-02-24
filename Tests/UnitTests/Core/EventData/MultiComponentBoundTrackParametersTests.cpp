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
#include "Acts/EventData/GenericBoundTrackParameters.hpp"
#include "Acts/EventData/MultiComponentTrackParameters.hpp"
#include "Acts/EventData/ParticleHypothesis.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Surfaces/CurvilinearSurface.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/TrackFitting/GsfOptions.hpp"

#include <initializer_list>
#include <memory>
#include <optional>
#include <tuple>
#include <vector>

using namespace Acts;

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

  const auto ap =
      MultiComponentBoundTrackParameters(surface, a, particleHypothesis);
  const auto bp =
      MultiComponentBoundTrackParameters(surface, b, particleHypothesis);
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
      std::vector<std::tuple<double, BoundVector, std::optional<BoundMatrix>>>
          a;
      for (int i = 0; i < 4; ++i) {
        a.emplace_back(0.25, single_pars.parameters(),
                       single_pars.covariance());
      }
      return MultiComponentBoundTrackParameters(surface, a, particleHypothesis);
    }();
    const auto multi_pars_as_single =
        multi_pars.merge(ComponentMergeMethod::eMean);

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
    if (cov && *cov != BoundMatrix::Zero()) {
      BOOST_CHECK_EQUAL(*multi_pars_as_single.covariance(),
                        *single_pars.covariance());
    } else {
      BOOST_CHECK(!multi_pars_as_single.covariance());
    }
  }
}
BOOST_AUTO_TEST_SUITE_END()

}  // namespace ActsTests

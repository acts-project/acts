// This file is part of the Acts project.
//
// Copyright (C) 2017-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include <Acts/EventData/Charge.hpp>
#include <Acts/EventData/MultiComponentBoundTrackParameters.hpp>
#include <Acts/Surfaces/PlaneSurface.hpp>
#include <Acts/Utilities/TrackParameterHelpers.hpp>

using namespace Acts;

BOOST_AUTO_TEST_CASE(test_constructors) {
  std::vector<std::tuple<double, BoundVector, BoundSymMatrix>> a;
  a.push_back({1.0, BoundVector::Ones(), BoundSymMatrix::Identity()});

  std::vector<std::tuple<double, BoundVector, std::optional<BoundSymMatrix>>> b;
  b.push_back({1.0, BoundVector::Ones(), BoundSymMatrix::Identity()});

  auto surface = Acts::Surface::makeShared<Acts::PlaneSurface>(
      Vector3::Ones(), Vector3::Ones().normalized());

  const auto ap = MultiComponentBoundTrackParameters(surface, a);
  const auto bp = MultiComponentBoundTrackParameters(surface, b);
  const auto aps = MultiComponentBoundTrackParameters(
      surface, std::get<1>(a.front()), std::get<2>(a.front()));
  const auto bps = MultiComponentBoundTrackParameters(
      surface, std::get<1>(b.front()), std::get<2>(b.front()));

  BOOST_CHECK(b == ap.components());
  BOOST_CHECK(ap.components() == bp.components());
  BOOST_CHECK(bp.components() == aps.components());
  BOOST_CHECK(aps.components() == bps.components());
}

BOOST_AUTO_TEST_CASE(test_accessors) {
  using cov_t = std::optional<BoundSymMatrix>;
  for (const auto &cov : {cov_t{}, cov_t{BoundSymMatrix::Identity()},
                          cov_t{BoundSymMatrix::Identity()}}) {
    auto surface = Acts::Surface::makeShared<Acts::PlaneSurface>(
        Vector3::Ones(), Vector3::Ones().normalized());

    const BoundTrackParameters single_pars(surface, BoundVector::Ones(), cov);

    const auto multi_pars = [&]() {
      std::vector<
          std::tuple<double, BoundVector, std::optional<BoundSymMatrix>>>
          a;
      for (int i = 0; i < 4; ++i) {
        a.push_back({0.25, single_pars.parameters(), single_pars.covariance()});
      }
      return MultiComponentBoundTrackParameters(surface, a);
    }();

    BOOST_CHECK_EQUAL(TrackParameterHelpers::absoluteMomentum(multi_pars, 1),
                      TrackParameterHelpers::absoluteMomentum(single_pars, 1));
    BOOST_CHECK_EQUAL(TrackParameterHelpers::charge(multi_pars, 1),
                      TrackParameterHelpers::charge(single_pars, 1));
    BOOST_CHECK_EQUAL(multi_pars.fourPosition(GeometryContext{}),
                      single_pars.fourPosition(GeometryContext{}));
    BOOST_CHECK_EQUAL(TrackParameterHelpers::momentum(multi_pars, 1),
                      TrackParameterHelpers::momentum(single_pars, 1));
    BOOST_CHECK_EQUAL(multi_pars.parameters(), single_pars.parameters());
    BOOST_CHECK_EQUAL(multi_pars.position(GeometryContext{}),
                      single_pars.position(GeometryContext{}));
    BOOST_CHECK_EQUAL(
        TrackParameterHelpers::transverseMomentum(multi_pars, 1),
        TrackParameterHelpers::transverseMomentum(single_pars, 1));
    BOOST_CHECK_EQUAL(multi_pars.direction(), single_pars.direction());

    // Check the behaviour for std::nullopt or zero covariance
    if (cov && *cov != BoundSymMatrix::Zero()) {
      BOOST_CHECK_EQUAL(*multi_pars.covariance(), *single_pars.covariance());
    } else {
      BOOST_CHECK(not multi_pars.covariance());
    }
  }
}

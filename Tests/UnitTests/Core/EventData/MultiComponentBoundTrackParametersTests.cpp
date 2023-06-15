// This file is part of the Acts project.
//
// Copyright (C) 2017-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/Definitions/Units.hpp"
#include "Acts/EventData/SingleBoundTrackParameters.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include <Acts/EventData/Charge.hpp>
#include <Acts/EventData/MultiComponentBoundTrackParameters.hpp>
#include <Acts/Surfaces/PlaneSurface.hpp>

#include <algorithm>
#include <initializer_list>
#include <memory>
#include <optional>
#include <tuple>
#include <vector>

using namespace Acts;

BOOST_AUTO_TEST_CASE(test_constructors) {
  std::vector<std::tuple<double, BoundVector, BoundSymMatrix>> a;
  a.push_back({1.0, BoundVector::Ones(), BoundSymMatrix::Identity()});

  std::vector<std::tuple<double, BoundVector, std::optional<BoundSymMatrix>>> b;
  b.push_back({1.0, BoundVector::Ones(), BoundSymMatrix::Identity()});

  auto surface = Acts::Surface::makeShared<Acts::PlaneSurface>(
      Vector3::Ones(), Vector3::Ones().normalized());

  const auto ap = MultiComponentBoundTrackParameters<SinglyCharged>(surface, a);
  const auto bp = MultiComponentBoundTrackParameters<SinglyCharged>(surface, b);
  const auto aps = MultiComponentBoundTrackParameters<SinglyCharged>(
      surface, std::get<1>(a.front()), std::get<2>(a.front()));
  const auto bps = MultiComponentBoundTrackParameters<SinglyCharged>(
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

    const SingleBoundTrackParameters<SinglyCharged> single_pars(
        surface, BoundVector::Ones(), cov);

    const auto multi_pars = [&]() {
      std::vector<
          std::tuple<double, BoundVector, std::optional<BoundSymMatrix>>>
          a;
      for (int i = 0; i < 4; ++i) {
        a.push_back({0.25, single_pars.parameters(), single_pars.covariance()});
      }
      return MultiComponentBoundTrackParameters<SinglyCharged>(surface, a);
    }();

    BOOST_CHECK_EQUAL(multi_pars.absoluteMomentum(),
                      single_pars.absoluteMomentum());
    BOOST_CHECK_EQUAL(multi_pars.charge(), single_pars.charge());
    BOOST_CHECK_EQUAL(multi_pars.fourPosition(GeometryContext{}),
                      single_pars.fourPosition(GeometryContext{}));
    BOOST_CHECK_EQUAL(multi_pars.momentum(), single_pars.momentum());
    BOOST_CHECK_EQUAL(multi_pars.parameters(), single_pars.parameters());
    BOOST_CHECK_EQUAL(multi_pars.position(GeometryContext{}),
                      single_pars.position(GeometryContext{}));
    BOOST_CHECK_EQUAL(multi_pars.transverseMomentum(),
                      single_pars.transverseMomentum());
    BOOST_CHECK_EQUAL(multi_pars.unitDirection(), single_pars.unitDirection());

    // Check the behaviour for std::nullopt or zero covariance
    if (cov && *cov != BoundSymMatrix::Zero()) {
      BOOST_CHECK_EQUAL(*multi_pars.covariance(), *single_pars.covariance());
    } else {
      BOOST_CHECK(not multi_pars.covariance());
    }
  }
}

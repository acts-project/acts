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

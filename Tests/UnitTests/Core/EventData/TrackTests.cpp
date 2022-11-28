// This file is part of the Acts project.
//
// Copyright (C) 2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/test/data/test_case.hpp>
#include <boost/test/unit_test.hpp>

#include "Acts/EventData/MultiTrajectory.hpp"
#include "Acts/EventData/Track.hpp"
#include "Acts/EventData/VectorMultiTrajectory.hpp"
#include "Acts/EventData/VectorTrackContainer.hpp"

namespace {

using namespace Acts;
using namespace Acts::HashedStringLiteral;
// using namespace Acts::UnitLiterals;
// using namespace Acts::Test;
using MultiTrajectoryTraits::IndexType;
namespace bd = boost::unit_test::data;
}  // namespace

BOOST_AUTO_TEST_SUITE(EventDataTrack)

BOOST_AUTO_TEST_CASE(Build) {
  VectorMultiTrajectory mtj{};
  TrackContainer tc{VectorTrackContainer{}, mtj};
  auto idx = tc.addTrack();
  auto t = tc.getTrack(idx);
  t.component<IndexType, "tipIndex"_hash>() = 5;

  BOOST_CHECK_EQUAL((t.component<IndexType, "tipIndex"_hash>()), 5);

  BoundVector pars;
  pars.setRandom();
  t.parameters() = pars;
  BOOST_CHECK_EQUAL(t.parameters(), pars);

  BoundMatrix cov;
  cov.setRandom();
  t.covariance() = cov;
  BOOST_CHECK_EQUAL(t.covariance(), cov);

  // const checks: should not compile
  // const auto& ctc = tc;
  // ctc.getTrack(idx).covariance().setRandom();
  // const auto& ctp = t;
  // ctp.covariance().setRandom();
}

BOOST_AUTO_TEST_CASE(BuildReadOnly) {
  ConstVectorMultiTrajectory mtj{};
  TrackContainer tc{ConstVectorTrackContainer{}, mtj};
}

BOOST_AUTO_TEST_CASE(DynamicColumns) {
  VectorMultiTrajectory mtj{};
  TrackContainer tc{VectorTrackContainer{}, mtj};
  BOOST_CHECK(!tc.hasColumn("col_a"_hash));
  tc.addColumn<float>("col_a");
  BOOST_CHECK(tc.hasColumn("col_a"_hash));

  auto t = tc.getTrack(tc.addTrack());
  t.component<float>("col_a") = 5.6f;
  BOOST_CHECK_EQUAL((t.component<float, "col_a"_hash>()), 5.6f);
}

BOOST_AUTO_TEST_SUITE_END()

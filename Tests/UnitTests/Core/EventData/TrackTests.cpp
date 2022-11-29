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
#include "Acts/EventData/TrackStatePropMask.hpp"
#include "Acts/EventData/VectorMultiTrajectory.hpp"
#include "Acts/EventData/VectorTrackContainer.hpp"
#include "Acts/Tests/CommonHelpers/TestTrackState.hpp"

namespace {

using namespace Acts;
using namespace Acts::HashedStringLiteral;
using namespace Acts::Test;
using MultiTrajectoryTraits::IndexType;
namespace bd = boost::unit_test::data;

const GeometryContext gctx;
// fixed seed for reproducible tests
std::default_random_engine rng(31415);
}  // namespace

BOOST_AUTO_TEST_SUITE(EventDataTrack)

BOOST_AUTO_TEST_CASE(Build) {
  VectorMultiTrajectory mtj{};
  VectorTrackContainer vtc{};
  TrackContainer tc{vtc, mtj};

  static_assert(
      std::is_same_v<decltype(tc), TrackContainer<VectorTrackContainer,
                                                  VectorMultiTrajectory>>,
      "Incorrect deduction");

  static_assert(!tc.ReadOnly, "Should be read only");

  auto idx = tc.addTrack();
  auto t = tc.getTrack(idx);
  t.component<IndexType, "tipIndex"_hash>() = 5;

  BOOST_CHECK_EQUAL((t.component<IndexType, "tipIndex"_hash>()), 5);
  BOOST_CHECK_EQUAL(t.tipIndex(), 5);
  t.tipIndex() = 6;
  BOOST_CHECK_EQUAL(t.tipIndex(), 6);

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

BOOST_AUTO_TEST_CASE(TrackStateAccess) {
  VectorMultiTrajectory traj;

  auto mkts = [&](auto prev) {
    if constexpr (std::is_same_v<decltype(prev), IndexType>) {
      auto ts =
          traj.getTrackState(traj.addTrackState(TrackStatePropMask::All, prev));
      TestTrackState pc(rng, 2u);
      fillTrackState(pc, TrackStatePropMask::All, ts);
      return ts;
    } else {
      auto ts = traj.getTrackState(
          traj.addTrackState(TrackStatePropMask::All, prev.index()));
      TestTrackState pc(rng, 2u);
      fillTrackState(pc, TrackStatePropMask::All, ts);
      return ts;
    }
  };

  auto ts1 = mkts(MultiTrajectoryTraits::kInvalid);
  auto ts2 = mkts(ts1);
  auto ts3 = mkts(ts2);
  auto ts4 = mkts(ts3);
  auto ts5 = mkts(ts4);

  VectorMultiTrajectory mtj{};
  VectorTrackContainer vtc{};
  TrackContainer tc{vtc, mtj};

  auto t = tc.getTrack(tc.addTrack());
  t.tipIndex() = ts5.index();

  t.visitTrackStates([](const auto& state) {});
}

BOOST_AUTO_TEST_CASE(BuildReadOnly) {
  ConstVectorMultiTrajectory mtj{};
  ConstVectorTrackContainer vtc{};
  TrackContainer tc{vtc, mtj};

  static_assert(
      std::is_same_v<decltype(tc), TrackContainer<ConstVectorTrackContainer,
                                                  ConstVectorMultiTrajectory>>,
      "Incorrect deduction");

  static_assert(tc.ReadOnly, "Should be read only");
}

BOOST_AUTO_TEST_CASE(DynamicColumns) {
  VectorMultiTrajectory mtj{};
  VectorTrackContainer vtc{};
  TrackContainer tc{vtc, mtj};
  BOOST_CHECK(!tc.hasColumn("col_a"_hash));
  tc.addColumn<float>("col_a");
  BOOST_CHECK(tc.hasColumn("col_a"_hash));

  auto t = tc.getTrack(tc.addTrack());
  t.component<float>("col_a") = 5.6f;
  BOOST_CHECK_EQUAL((t.component<float, "col_a"_hash>()), 5.6f);
}

BOOST_AUTO_TEST_SUITE_END()

// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
#include <boost/test/unit_test.hpp>

#include "Acts/EventData/VectorMultiTrajectory.hpp"
#include "Acts/EventData/VectorTrackContainer.hpp"
#include "Acts/Utilities/Zip.hpp"

#include <numeric>

using namespace Acts;
using namespace Acts::HashedStringLiteral;
using MultiTrajectoryTraits::IndexType;

BOOST_AUTO_TEST_SUITE(EventDataTrack)

BOOST_AUTO_TEST_CASE(CopyTracksIncludingDynamicColumns) {
  // mutable source
  VectorTrackContainer vtc{};
  VectorMultiTrajectory mtj{};
  TrackContainer tc{vtc, mtj};
  tc.addColumn<size_t>("counter");
  tc.addColumn<bool>("odd");

  TrackContainer tc2{VectorTrackContainer{}, VectorMultiTrajectory{}};
  // doesn't have the dynamic column

  TrackContainer tc3{VectorTrackContainer{}, VectorMultiTrajectory{}};
  tc3.addColumn<size_t>("counter");
  tc3.addColumn<bool>("odd");

  for (size_t i = 0; i < 10; i++) {
    auto t = tc.getTrack(tc.addTrack());
    auto ts = t.appendTrackState();
    ts.predicted() = BoundVector::Ones();
    ts = t.appendTrackState();
    ts.predicted().setOnes();
    ts.predicted() *= 2;

    ts = t.appendTrackState();
    ts.predicted().setOnes();
    ts.predicted() *= 3;

    t.template component<size_t>("counter") = i;
    t.template component<bool>("odd") = i % 2 == 0;

    auto t2 = tc2.getTrack(tc2.addTrack());
    BOOST_CHECK_THROW(t2.copyFrom(t),
                      std::invalid_argument);  // this should fail

    auto t3 = tc3.getTrack(tc3.addTrack());
    t3.copyFrom(t);  // this should work

    BOOST_CHECK_NE(t3.tipIndex(), MultiTrajectoryTraits::kInvalid);
    BOOST_CHECK(t3.nTrackStates() > 0);
    BOOST_REQUIRE_EQUAL(t.nTrackStates(), t3.nTrackStates());

    for (auto [tsa, tsb] :
         zip(t.trackStatesReversed(), t3.trackStatesReversed())) {
      BOOST_CHECK_EQUAL(tsa.predicted(), tsb.predicted());
    }

    BOOST_CHECK_EQUAL(t.template component<size_t>("counter"),
                      t3.template component<size_t>("counter"));
    BOOST_CHECK_EQUAL(t.template component<bool>("odd"),
                      t3.template component<bool>("odd"));
  }

  size_t before = mtj.size();
  TrackContainer tc4{ConstVectorTrackContainer{vtc},
                     ConstVectorMultiTrajectory{mtj}};

  BOOST_REQUIRE_EQUAL(tc4.trackStateContainer().size(), before);

  TrackContainer tc5{VectorTrackContainer{}, VectorMultiTrajectory{}};
  tc5.addColumn<size_t>("counter");
  tc5.addColumn<bool>("odd");

  for (size_t i = 0; i < 10; i++) {
    auto t4 = tc4.getTrack(i);  // const source!
    BOOST_CHECK_NE(t4.nTrackStates(), 0);

    auto t5 = tc5.getTrack(tc5.addTrack());
    t5.copyFrom(t4);  // this should work

    BOOST_CHECK_NE(t5.tipIndex(), MultiTrajectoryTraits::kInvalid);
    BOOST_CHECK(t5.nTrackStates() > 0);
    BOOST_REQUIRE_EQUAL(t4.nTrackStates(), t5.nTrackStates());

    for (auto [tsa, tsb] :
         zip(t4.trackStatesReversed(), t5.trackStatesReversed())) {
      BOOST_CHECK_EQUAL(tsa.predicted(), tsb.predicted());
    }

    BOOST_CHECK_EQUAL(t4.template component<size_t>("counter"),
                      t5.template component<size_t>("counter"));
    BOOST_CHECK_EQUAL(t4.template component<bool>("odd"),
                      t5.template component<bool>("odd"));
  }
}

BOOST_AUTO_TEST_CASE(ReverseTrackStates) {
  VectorTrackContainer vtc{};
  VectorMultiTrajectory mtj{};
  TrackContainer tc{vtc, mtj};

  auto t = tc.getTrack(tc.addTrack());

  for (size_t i = 0; i < 4; i++) {
    auto ts = t.appendTrackState();
    ts.jacobian() = Acts::BoundMatrix::Identity() * i;
  }

  std::vector<IndexType> exp;
  exp.resize(t.nTrackStates());
  std::iota(exp.rbegin(), exp.rend(), 0);
  std::vector<IndexType> act;
  std::transform(t.trackStatesReversed().begin(), t.trackStatesReversed().end(),
                 std::back_inserter(act),
                 [](const auto& ts) { return ts.index(); });

  // jacobians count up
  for (const auto [e, ts] : zip(exp, t.trackStatesReversed())) {
    BOOST_CHECK_EQUAL(ts.jacobian(), Acts::BoundMatrix::Identity() * e);
  }

  BOOST_CHECK_EQUAL_COLLECTIONS(exp.begin(), exp.end(), act.begin(), act.end());

  // reverse!
  t.reverseTrackStates();

  std::iota(exp.begin(), exp.end(), 0);
  act.clear();
  std::transform(t.trackStatesReversed().begin(), t.trackStatesReversed().end(),
                 std::back_inserter(act),
                 [](const auto& ts) { return ts.index(); });
  BOOST_CHECK_EQUAL_COLLECTIONS(exp.begin(), exp.end(), act.begin(), act.end());

  // jacobians stay with their track states
  for (const auto [e, ts] : zip(exp, t.trackStatesReversed())) {
    BOOST_CHECK_EQUAL(ts.jacobian(), Acts::BoundMatrix::Identity() * e);
  }

  // back to original!
  t.reverseTrackStates();

  // jacobians stay with their track states
  for (const auto [e, ts] : zip(exp, t.trackStates())) {
    BOOST_CHECK_EQUAL(ts.jacobian(), Acts::BoundMatrix::Identity() * e);
  }

  // reverse with jacobians
  t.reverseTrackStates(true);

  std::reverse(exp.begin(), exp.end());
  std::rotate(exp.rbegin(), std::next(exp.rbegin()), exp.rend());

  for (const auto [e, ts] : zip(exp, t.trackStates())) {
    Acts::BoundMatrix expJac;
    if (e == 0) {
      expJac = Acts::BoundMatrix::Zero();
    } else {
      expJac = (Acts::BoundMatrix::Identity() * e).inverse();
    }

    BOOST_CHECK_EQUAL(ts.jacobian(), expJac);
  }

  // now back to original order, revert jacobians again
  t.reverseTrackStates(true);

  // reset exp to range(0, N)
  std::iota(exp.begin(), exp.end(), 0);

  for (const auto [e, ts] : zip(exp, t.trackStates())) {
    BOOST_CHECK_EQUAL(ts.jacobian(), Acts::BoundMatrix::Identity() * e);
  }
}

BOOST_AUTO_TEST_SUITE_END()

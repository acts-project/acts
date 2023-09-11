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

using namespace Acts;

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
BOOST_AUTO_TEST_SUITE_END()

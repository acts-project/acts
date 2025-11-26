// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/EventData/AnyTrackState.hpp"
#include "Acts/EventData/MultiTrajectory.hpp"
#include "Acts/EventData/VectorMultiTrajectory.hpp"
#include "Acts/EventData/VectorTrackContainer.hpp"
#include "Acts/EventData/TrackContainer.hpp"
#include "Acts/EventData/TrackStateProxy.hpp"

namespace {

using namespace Acts;

struct TestTrackStateFixture {
  using Trajectory = VectorMultiTrajectory;
  using TrackContainerBackend = VectorTrackContainer;
  using Container = TrackContainer<TrackContainerBackend, Trajectory,
                                   detail::RefHolder>;

  TestTrackStateFixture() : container(trackContainer, trajectory) {}

  Container container;
  TrackContainerBackend trackContainer;
  Trajectory trajectory;
};

}  // namespace

BOOST_AUTO_TEST_SUITE(EventDataAnyTrackState)

BOOST_FIXTURE_TEST_CASE(ConstructDefault, TestTrackStateFixture) {
  AnyMutableTrackState ts;
  BOOST_CHECK(!ts.isValid());
  BOOST_CHECK(!ts);
}

BOOST_FIXTURE_TEST_CASE(WrapTrackStateProxy, TestTrackStateFixture) {
  auto track = container.makeTrack();
  auto state = container.trackStateContainer().makeTrackState();
  track.tipIndex() = state.index();

  auto proxy = container.trackStateContainer().getTrackState(state.index());
  AnyMutableTrackState anyState(proxy);

  BOOST_CHECK(anyState.isValid());
  BOOST_CHECK_EQUAL(anyState.index(), proxy.index());
}

BOOST_FIXTURE_TEST_CASE(AccessFiltered, TestTrackStateFixture) {
  auto track = container.makeTrack();
  auto state = container.trackStateContainer().makeTrackState();
  track.tipIndex() = state.index();

  state.predicted() = ActsVector<eBoundSize>::Ones();
  state.filtered() = 2. * ActsVector<eBoundSize>::Ones();

  AnyMutableTrackState anyState(state);

  BOOST_CHECK_EQUAL(anyState.predicted()[eBoundLoc0], 1.);
  BOOST_CHECK_EQUAL(anyState.filtered()[eBoundLoc0], 2.);

  auto filtered = anyState.filtered();
  filtered[eBoundLoc0] = 3.;
  BOOST_CHECK_EQUAL(state.filtered()[eBoundLoc0], 3.);
}

BOOST_FIXTURE_TEST_CASE(AccessCalibratedFixedSize, TestTrackStateFixture) {
  auto track = container.makeTrack();
  auto state = container.trackStateContainer().makeTrackState();
  track.tipIndex() = state.index();

  state.allocateCalibrated(2);
  auto calibrated = state.template calibrated<2>();
  calibrated[0] = 1.5;
  calibrated[1] = 2.5;
  auto calibratedCov = state.template calibratedCovariance<2>();
  calibratedCov.setZero();
  calibratedCov(0, 0) = 0.1;
  calibratedCov(1, 1) = 0.2;

  AnyMutableTrackState anyState(state);

  auto view = anyState.calibrated<2>();
  BOOST_CHECK_CLOSE(view[0], 1.5, 1e-6);
  BOOST_CHECK_CLOSE(view[1], 2.5, 1e-6);

  view[0] = 4.5;
  BOOST_CHECK_CLOSE(state.calibrated<2>()[0], 4.5, 1e-6);

  auto covView = anyState.calibratedCovariance<2>();
  BOOST_CHECK_CLOSE(covView(0, 0), 0.1, 1e-6);
  BOOST_CHECK_CLOSE(covView(1, 1), 0.2, 1e-6);
  covView(0, 0) = 0.5;
  BOOST_CHECK_CLOSE(state.calibratedCovariance<2>()(0, 0), 0.5, 1e-6);

  AnyConstTrackState constState(state);
  auto constView = constState.calibrated<2>();
  BOOST_CHECK_CLOSE(constView[0], 4.5, 1e-6);
  auto constCov = constState.calibratedCovariance<2>();
  BOOST_CHECK_CLOSE(constCov(0, 0), 0.5, 1e-6);
}

BOOST_AUTO_TEST_SUITE_END()

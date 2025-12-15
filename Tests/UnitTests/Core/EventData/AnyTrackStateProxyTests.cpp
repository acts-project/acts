// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/EventData/AnyTrackStateProxy.hpp"
#include "Acts/EventData/ProxyAccessor.hpp"
#include "Acts/EventData/TrackContainer.hpp"
#include "Acts/EventData/TrackStateProxy.hpp"
#include "Acts/EventData/VectorMultiTrajectory.hpp"
#include "Acts/EventData/VectorTrackContainer.hpp"
#include "Acts/Utilities/HashedString.hpp"

namespace {

using namespace Acts;
using namespace Acts::HashedStringLiteral;

struct TestTrackStateFixture {
  using Trajectory = VectorMultiTrajectory;
  using TrackContainerBackend = VectorTrackContainer;
  using Container =
      TrackContainer<TrackContainerBackend, Trajectory, detail::RefHolder>;

  TestTrackStateFixture() : container(trackContainer, trajectory) {}

  Container container;
  TrackContainerBackend trackContainer;
  Trajectory trajectory;
};

}  // namespace

BOOST_AUTO_TEST_SUITE(EventDataAnyTrackState)

BOOST_FIXTURE_TEST_CASE(WrapTrackStateProxy, TestTrackStateFixture) {
  auto track = container.makeTrack();
  auto state = container.trackStateContainer().makeTrackState();
  track.tipIndex() = state.index();

  auto proxy = container.trackStateContainer().getTrackState(state.index());
  AnyMutableTrackStateProxy anyState(proxy);

  BOOST_CHECK_EQUAL(anyState.index(), proxy.index());
}

BOOST_FIXTURE_TEST_CASE(AccessFiltered, TestTrackStateFixture) {
  auto track = container.makeTrack();
  auto state = container.trackStateContainer().makeTrackState();
  track.tipIndex() = state.index();

  state.predicted() = ActsVector<eBoundSize>::Ones();
  state.filtered() = 2. * ActsVector<eBoundSize>::Ones();

  AnyMutableTrackStateProxy anyState(state);

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

  AnyMutableTrackStateProxy anyState(state);

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

  AnyConstTrackStateProxy constState(state);
  auto constView = constState.calibrated<2>();
  BOOST_CHECK_CLOSE(constView[0], 4.5, 1e-6);
  auto constCov = constState.calibratedCovariance<2>();
  BOOST_CHECK_CLOSE(constCov(0, 0), 0.5, 1e-6);
}

BOOST_FIXTURE_TEST_CASE(AccessEffectiveCalibratedDynamic,
                        TestTrackStateFixture) {
  auto track = container.makeTrack();
  auto state = container.trackStateContainer().makeTrackState();
  track.tipIndex() = state.index();

  constexpr std::size_t measdim = 3u;
  state.allocateCalibrated(measdim);
  {
    auto dyn = state.effectiveCalibrated();
    for (std::size_t i = 0; i < measdim; ++i) {
      dyn[i] = 1. + static_cast<double>(i);
    }
  }
  {
    auto dynCov = state.effectiveCalibratedCovariance();
    dynCov.setZero();
    for (std::size_t i = 0; i < measdim; ++i) {
      dynCov(i, i) = 0.1 * static_cast<double>(i + 1);
    }
  }

  AnyMutableTrackStateProxy anyState(state);
  BOOST_CHECK_EQUAL(anyState.calibratedSize(), measdim);
  auto eff = anyState.effectiveCalibrated();
  BOOST_CHECK_EQUAL(eff.size(), measdim);
  BOOST_CHECK_CLOSE(eff[0], 1., 1e-6);
  BOOST_CHECK_CLOSE(eff[2], 3., 1e-6);
  eff[1] = 5.5;
  BOOST_CHECK_CLOSE(state.effectiveCalibrated()[1], 5.5, 1e-6);

  auto effCov = anyState.effectiveCalibratedCovariance();
  BOOST_CHECK_EQUAL(effCov.rows(), measdim);
  BOOST_CHECK_EQUAL(effCov.cols(), measdim);
  BOOST_CHECK_CLOSE(effCov(0, 0), 0.1, 1e-6);
  BOOST_CHECK_CLOSE(effCov(2, 2), 0.3, 1e-6);
  effCov(1, 1) = 1.7;
  BOOST_CHECK_CLOSE(state.effectiveCalibratedCovariance()(1, 1), 1.7, 1e-6);

  AnyConstTrackStateProxy constState(state);
  BOOST_CHECK_EQUAL(constState.calibratedSize(), measdim);
  auto constEff = constState.effectiveCalibrated();
  BOOST_CHECK_EQUAL(constEff.size(), measdim);
  BOOST_CHECK_CLOSE(constEff[1], 5.5, 1e-6);
  auto constEffCov = constState.effectiveCalibratedCovariance();
  BOOST_CHECK_EQUAL(constEffCov.rows(), measdim);
  BOOST_CHECK_CLOSE(constEffCov(1, 1), 1.7, 1e-6);
}

BOOST_FIXTURE_TEST_CASE(ProxyAccessorWithAnyTrackState, TestTrackStateFixture) {
  container.trackStateContainer().addColumn<float>("customFloat");

  auto track = container.makeTrack();
  auto state = container.trackStateContainer().makeTrackState();
  track.tipIndex() = state.index();

  state.template component<float>("customFloat"_hash) = 0.25f;

  ProxyAccessor<float> mutableAccessor("customFloat");
  ConstProxyAccessor<float> constAccessor("customFloat");

  AnyMutableTrackStateProxy anyState(state);
  BOOST_CHECK_CLOSE(mutableAccessor(anyState), 0.25f, 1e-6);
  mutableAccessor(anyState) = 0.75f;
  BOOST_CHECK_CLOSE(state.template component<float>("customFloat"_hash), 0.75f,
                    1e-6);

  AnyConstTrackStateProxy constState(state);
  BOOST_CHECK_CLOSE(constAccessor(constState), 0.75f, 1e-6);
  BOOST_CHECK(constAccessor.hasColumn(constState));
}

BOOST_AUTO_TEST_SUITE_END()

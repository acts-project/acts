// This file is part of the ACTS project.
//
// Copyright (C) 2025 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/EventData/MultiTrajectory.hpp"
#include "Acts/EventData/TrackContainer.hpp"
#include "Acts/EventData/TrackHandle.hpp"
#include "Acts/EventData/TrackProxy.hpp"
#include "Acts/EventData/VectorMultiTrajectory.hpp"
#include "Acts/EventData/VectorTrackContainer.hpp"
#include "Acts/EventData/detail/TestTrackState.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Surfaces/PerigeeSurface.hpp"
#include "Acts/Utilities/HashedString.hpp"
#include "Acts/Utilities/Holders.hpp"

#include <memory>
#include <random>
#include <vector>

namespace {

using namespace Acts;
using namespace Acts::HashedStringLiteral;
using namespace Acts::detail::Test;

const GeometryContext gctx;
std::default_random_engine rng(31415);

}  // namespace

namespace ActsTests {

BOOST_AUTO_TEST_SUITE(EventDataSuite)

BOOST_AUTO_TEST_CASE(TrackHandleConstruction) {
  // Create a track container
  VectorTrackContainer vtc;
  VectorMultiTrajectory mtj;
  TrackContainer tc{vtc, mtj};

  auto t = tc.makeTrack();
  t.tipIndex() = 42;
  t.stemIndex() = 24;

  // Create a handle from the track
  TrackHandle handle(t);

  BOOST_CHECK(handle.isValid());
  BOOST_CHECK(static_cast<bool>(handle));
  BOOST_CHECK_EQUAL(handle.tipIndex(), 42);
  BOOST_CHECK_EQUAL(handle.stemIndex(), 24);
}

BOOST_AUTO_TEST_CASE(TrackHandleConstConstruction) {
  VectorTrackContainer vtc;
  VectorMultiTrajectory mtj;
  TrackContainer tc{vtc, mtj};

  auto t = tc.makeTrack();
  t.tipIndex() = 100;

  const auto& ct = t;
  TrackHandle handle(ct);

  BOOST_CHECK(handle.isValid());
  BOOST_CHECK_EQUAL(handle.tipIndex(), 100);
}

BOOST_AUTO_TEST_CASE(TrackHandleEmptyConstruction) {
  TrackHandle handle;

  BOOST_CHECK(!handle.isValid());
  BOOST_CHECK(!static_cast<bool>(handle));
}

BOOST_AUTO_TEST_CASE(TrackHandleBasicInterface) {
  VectorTrackContainer vtc;
  VectorMultiTrajectory mtj;
  TrackContainer tc{vtc, mtj};

  auto t = tc.makeTrack();
  t.tipIndex() = 10;
  t.stemIndex() = 5;

  // Set up track states
  auto ts1 = t.appendTrackState();
  ts1.predicted().setRandom();
  ts1.typeFlags().set(TrackStateFlag::MeasurementFlag);

  auto ts2 = t.appendTrackState();
  ts2.predicted().setRandom();
  ts2.typeFlags().set(TrackStateFlag::MeasurementFlag);

  auto ts3 = t.appendTrackState();
  ts3.predicted().setRandom();
  ts3.typeFlags().set(TrackStateFlag::OutlierFlag);

  TrackHandle handle(t);

  BOOST_CHECK_EQUAL(handle.tipIndex(), 10);
  BOOST_CHECK_EQUAL(handle.stemIndex(), 5);
  BOOST_CHECK_EQUAL(handle.nTrackStates(), 3);
  BOOST_CHECK_EQUAL(handle.nMeasurements(), 2);
  BOOST_CHECK_EQUAL(handle.nOutliers(), 1);
  BOOST_CHECK_EQUAL(handle.nHoles(), 0);
}

BOOST_AUTO_TEST_CASE(TrackHandleParameters) {
  VectorTrackContainer vtc;
  VectorMultiTrajectory mtj;
  TrackContainer tc{vtc, mtj};

  auto t = tc.makeTrack();

  auto perigee = Surface::makeShared<PerigeeSurface>(Vector3::Zero());
  t.setReferenceSurface(perigee);

  BoundVector params;
  params << 1.0, 2.0, 0.1, M_PI / 4.0, 0.5, 0.0;
  t.parameters() = params;

  BoundSquareMatrix cov = BoundSquareMatrix::Identity();
  t.covariance() = cov;

  t.particleHypothesis() = ParticleHypothesis::pion();

  TrackHandle handle(t);

  BOOST_CHECK(handle.hasReferenceSurface());
  BOOST_CHECK_EQUAL(&handle.referenceSurface(), perigee.get());
  BOOST_CHECK_EQUAL(handle.parameters().parameters(), params);
  BOOST_CHECK_EQUAL(handle.covariance(), cov);
  BOOST_CHECK_EQUAL(handle.particleHypothesis(), ParticleHypothesis::pion());
}

BOOST_AUTO_TEST_CASE(TrackHandleStatistics) {
  VectorTrackContainer vtc;
  VectorMultiTrajectory mtj;
  TrackContainer tc{vtc, mtj};

  auto t = tc.makeTrack();

  // Create track states with different flags
  auto ts1 = t.appendTrackState();
  ts1.typeFlags().set(TrackStateFlag::MeasurementFlag);
  ts1.chi2() = 2.5;

  auto ts2 = t.appendTrackState();
  ts2.typeFlags().set(TrackStateFlag::MeasurementFlag);
  ts2.chi2() = 3.5;

  auto ts3 = t.appendTrackState();
  ts3.typeFlags().set(TrackStateFlag::HoleFlag);

  auto ts4 = t.appendTrackState();
  ts4.typeFlags().set(TrackStateFlag::OutlierFlag);

  // Set shared hit
  auto ts5 = t.appendTrackState();
  ts5.typeFlags().set(TrackStateFlag::MeasurementFlag);
  ts5.typeFlags().set(TrackStateFlag::SharedHitFlag);

  TrackHandle handle(t);

  BOOST_CHECK_EQUAL(handle.nTrackStates(), 5);
  BOOST_CHECK_EQUAL(handle.nMeasurements(), 3);
  BOOST_CHECK_EQUAL(handle.nHoles(), 1);
  BOOST_CHECK_EQUAL(handle.nOutliers(), 1);
  BOOST_CHECK_EQUAL(handle.nSharedHits(), 1);
  BOOST_CHECK_CLOSE(handle.chi2(), 6.0, 0.001);  // 2.5 + 3.5
}

BOOST_AUTO_TEST_CASE(TrackHandleDynamicColumns) {
  VectorTrackContainer vtc;
  VectorMultiTrajectory mtj;
  TrackContainer tc{vtc, mtj};

  tc.addColumn<float>("custom_score");

  auto t = tc.makeTrack();
  t.template component<float>("custom_score") = 3.14f;

  TrackHandle handle(t);

  auto scoreAny = handle.component("custom_score"_hash);
  BOOST_CHECK(scoreAny.has_value());

  float score = handle.componentAs<float>("custom_score"_hash);
  BOOST_CHECK_CLOSE(score, 3.14f, 0.001);
}

BOOST_AUTO_TEST_CASE(TrackHandleNoOwnership) {
  VectorTrackContainer vtc;
  VectorMultiTrajectory mtj;
  TrackContainer tc{vtc, mtj};

  auto t = tc.makeTrack();
  t.tipIndex() = 99;

  // Create handle in a scope
  TrackIndexType tipFromHandle;
  {
    TrackHandle handle(t);
    tipFromHandle = handle.tipIndex();
    BOOST_CHECK_EQUAL(tipFromHandle, 99);
  }  // handle goes out of scope

  // Track should still be valid
  BOOST_CHECK_EQUAL(t.tipIndex(), 99);
  BOOST_CHECK_EQUAL(tc.size(), 1);
}

BOOST_AUTO_TEST_CASE(TrackHandleMultipleTypes) {
  // Test with RefHolder
  {
    VectorTrackContainer vtc;
    VectorMultiTrajectory mtj;
    TrackContainer<VectorTrackContainer, VectorMultiTrajectory,
                   detail::RefHolder>
        tc{vtc, mtj};

    auto t = tc.makeTrack();
    t.tipIndex() = 1;

    TrackHandle handle(t);
    BOOST_CHECK_EQUAL(handle.tipIndex(), 1);
  }

  // Test with ValueHolder
  {
    TrackContainer<VectorTrackContainer, VectorMultiTrajectory,
                   detail::ValueHolder>
        tc{VectorTrackContainer{}, VectorMultiTrajectory{}};

    auto t = tc.makeTrack();
    t.tipIndex() = 2;

    TrackHandle handle(t);
    BOOST_CHECK_EQUAL(handle.tipIndex(), 2);
  }

  // Test with shared_ptr
  {
    auto vtc = std::make_shared<VectorTrackContainer>();
    auto mtj = std::make_shared<VectorMultiTrajectory>();
    TrackContainer<VectorTrackContainer, VectorMultiTrajectory, std::shared_ptr>
        tc{vtc, mtj};

    auto t = tc.makeTrack();
    t.tipIndex() = 3;

    TrackHandle handle(t);
    BOOST_CHECK_EQUAL(handle.tipIndex(), 3);
  }
}

BOOST_AUTO_TEST_CASE(TrackHandleConstTrackContainer) {
  VectorTrackContainer mutVtc;
  VectorMultiTrajectory mutMtj;

  {
    TrackContainer mutTc{mutVtc, mutMtj};
    auto t = mutTc.makeTrack();
    t.tipIndex() = 123;
    t.appendTrackState();
  }

  ConstVectorTrackContainer vtc{std::move(mutVtc)};
  ConstVectorMultiTrajectory mtj{std::move(mutMtj)};

  TrackContainer tc{vtc, mtj};
  auto t = tc.getTrack(0);

  TrackHandle handle(t);

  BOOST_CHECK(handle.isValid());
  BOOST_CHECK_EQUAL(handle.tipIndex(), 123);
  BOOST_CHECK_EQUAL(handle.nTrackStates(), 1);
}

BOOST_AUTO_TEST_CASE(TrackHandleTypeErasure) {
  VectorTrackContainer vtc1;
  VectorMultiTrajectory mtj1;
  TrackContainer tc1{vtc1, mtj1};

  auto t1 = tc1.makeTrack();
  t1.tipIndex() = 10;
  t1.appendTrackState();

  auto vtc2 = std::make_shared<VectorTrackContainer>();
  auto mtj2 = std::make_shared<VectorMultiTrajectory>();
  TrackContainer<VectorTrackContainer, VectorMultiTrajectory, std::shared_ptr>
      tc2{vtc2, mtj2};

  auto t2 = tc2.makeTrack();
  t2.tipIndex() = 20;
  t2.appendTrackState();
  t2.appendTrackState();

  // Different concrete types, but can both be stored in TrackHandle
  std::vector<TrackHandle> handles;
  handles.emplace_back(t1);
  handles.emplace_back(t2);

  BOOST_CHECK_EQUAL(handles.size(), 2);
  BOOST_CHECK_EQUAL(handles[0].tipIndex(), 10);
  BOOST_CHECK_EQUAL(handles[0].nTrackStates(), 1);
  BOOST_CHECK_EQUAL(handles[1].tipIndex(), 20);
  BOOST_CHECK_EQUAL(handles[1].nTrackStates(), 2);
}

BOOST_AUTO_TEST_CASE(TrackHandleSize) {
  // Verify that TrackHandle is small (two pointers)
  BOOST_CHECK_EQUAL(sizeof(TrackHandle), 2 * sizeof(void*));
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace ActsTests

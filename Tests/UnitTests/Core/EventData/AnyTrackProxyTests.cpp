// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/EventData/AnyTrackProxy.hpp"
#include "Acts/EventData/MultiTrajectory.hpp"
#include "Acts/EventData/ProxyAccessor.hpp"
#include "Acts/EventData/TrackContainer.hpp"
#include "Acts/EventData/TrackProxy.hpp"
#include "Acts/EventData/VectorMultiTrajectory.hpp"
#include "Acts/EventData/VectorTrackContainer.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Surfaces/PerigeeSurface.hpp"
#include "Acts/Utilities/HashedString.hpp"
#include "Acts/Utilities/Holders.hpp"

#include <memory>
#include <vector>

namespace {

using namespace Acts::UnitLiterals;
using namespace Acts;
using namespace Acts::HashedStringLiteral;

const GeometryContext gctx;

// Helper to create a test track with some data
template <typename track_container_t>
void fillTestTrack(typename track_container_t::TrackProxy track) {
  auto surface = Acts::Surface::makeShared<PerigeeSurface>(Vector3::Zero());

  // Set reference surface
  track.setReferenceSurface(surface);

  // Set parameters
  auto params = track.parameters();
  params[eBoundLoc0] = 1.0;
  params[eBoundLoc1] = 2.0;
  params[eBoundTime] = 3.0;
  params[eBoundPhi] = 0.5;
  params[eBoundTheta] = std::numbers::pi / 4.0;
  params[eBoundQOverP] = 0.1;

  // Set covariance
  auto cov = track.covariance();
  cov.setIdentity();
  cov(eBoundLoc0, eBoundLoc0) = 0.01;
  cov(eBoundQOverP, eBoundQOverP) = 0.0001;

  // Set statistics
  track.chi2() = 12.5f;
  track.nDoF() = 10u;
  track.nMeasurements() = 8u;
  track.nHoles() = 1u;
  track.nOutliers() = 0u;
  track.nSharedHits() = 2u;

  // Set particle hypothesis
  track.setParticleHypothesis(ParticleHypothesis::pion());
}

}  // namespace

BOOST_AUTO_TEST_SUITE(EventDataAnyTrack)

BOOST_AUTO_TEST_CASE(ConstructFromTrackProxy_ValueHolder) {
  TrackContainer tc{VectorTrackContainer{}, VectorMultiTrajectory{}};

  auto track = tc.makeTrack();
  fillTestTrack<decltype(tc)>(track);

  // Construct from track proxy
  AnyMutableTrackProxy anyTrack(track);

  // Verify index
  BOOST_CHECK_EQUAL(anyTrack.index(), track.index());
}

BOOST_AUTO_TEST_CASE(ConstructFromTrackProxy_SharedPtr) {
  auto vtc = std::make_shared<VectorTrackContainer>();
  auto mtj = std::make_shared<VectorMultiTrajectory>();
  TrackContainer tc{vtc, mtj};

  auto track = tc.makeTrack();
  fillTestTrack<decltype(tc)>(track);

  // Construct from track proxy
  AnyMutableTrackProxy anyTrack(track);

  // Verify index
  BOOST_CHECK_EQUAL(anyTrack.index(), track.index());
}

BOOST_AUTO_TEST_CASE(ConstructFromConstTrackProxy) {
  VectorTrackContainer vtc;
  VectorMultiTrajectory mtj;
  TrackContainer tc{vtc, mtj};

  auto track = tc.makeTrack();
  fillTestTrack<decltype(tc)>(track);

  // Get const track proxy
  auto constTrack = tc.getTrack(track.index());

  // Construct from const track proxy
  AnyConstTrackProxy anyTrack(constTrack);

  BOOST_CHECK_EQUAL(anyTrack.index(), track.index());
}

BOOST_AUTO_TEST_CASE(AccessIndices) {
  VectorTrackContainer vtc;
  VectorMultiTrajectory mtj;
  TrackContainer tc{vtc, mtj};

  auto track = tc.makeTrack();
  fillTestTrack<decltype(tc)>(track);

  track.tipIndex() = 42;
  track.stemIndex() = 7;

  AnyMutableTrackProxy anyTrack(track);

  BOOST_CHECK_EQUAL(anyTrack.tipIndex(), 42u);
  BOOST_CHECK_EQUAL(anyTrack.stemIndex(), 7u);
  BOOST_CHECK_EQUAL(anyTrack.index(), track.index());
}

BOOST_AUTO_TEST_CASE(AccessReferenceSurface) {
  VectorTrackContainer vtc;
  VectorMultiTrajectory mtj;
  TrackContainer tc{vtc, mtj};

  auto track = tc.makeTrack();
  auto surface = Acts::Surface::makeShared<PerigeeSurface>(Vector3::Zero());
  track.setReferenceSurface(surface);

  AnyMutableTrackProxy anyTrack(track);

  BOOST_CHECK(anyTrack.hasReferenceSurface());
  BOOST_CHECK_EQUAL(&anyTrack.referenceSurface(), surface.get());
}

BOOST_AUTO_TEST_CASE(AccessParameters) {
  VectorTrackContainer vtc;
  VectorMultiTrajectory mtj;
  TrackContainer tc{vtc, mtj};

  auto track = tc.makeTrack();
  fillTestTrack<decltype(tc)>(track);

  AnyMutableTrackProxy anyTrack(track);

  // Check individual parameters
  BOOST_CHECK_CLOSE(anyTrack.parameter(eBoundLoc0), 1.0, 1e-6);
  BOOST_CHECK_CLOSE(anyTrack.parameter(eBoundLoc1), 2.0, 1e-6);
  BOOST_CHECK_CLOSE(anyTrack.parameter(eBoundTime), 3.0, 1e-6);
  BOOST_CHECK_CLOSE(anyTrack.parameter(eBoundPhi), 0.5, 1e-6);
  BOOST_CHECK_CLOSE(anyTrack.parameter(eBoundTheta), std::numbers::pi / 4.0,
                    1e-6);
  BOOST_CHECK_CLOSE(anyTrack.parameter(eBoundQOverP), 0.1, 1e-6);

  // Check convenience accessors
  BOOST_CHECK_CLOSE(anyTrack.phi(), 0.5, 1e-6);
  BOOST_CHECK_CLOSE(anyTrack.theta(), std::numbers::pi / 4.0, 1e-6);
  BOOST_CHECK_CLOSE(anyTrack.qOverP(), 0.1, 1e-6);

  // Verify parameters() map provides reference semantics
  auto paramsView = anyTrack.parameters();
  paramsView[eBoundLoc0] = 9.0;
  BOOST_CHECK_CLOSE(track.parameters()[eBoundLoc0], 9.0, 1e-6);

  AnyConstTrackProxy constTrack(track);
  auto constParamsView = constTrack.parameters();
  BOOST_CHECK_CLOSE(constParamsView[eBoundLoc0], 9.0, 1e-6);
}

BOOST_AUTO_TEST_CASE(AccessCovariance) {
  VectorTrackContainer vtc;
  VectorMultiTrajectory mtj;
  TrackContainer tc{vtc, mtj};

  auto track = tc.makeTrack();
  fillTestTrack<decltype(tc)>(track);

  AnyMutableTrackProxy anyTrack(track);

  // Check covariance elements
  BOOST_CHECK_CLOSE(anyTrack.covariance(eBoundLoc0, eBoundLoc0), 0.01, 1e-6);
  BOOST_CHECK_CLOSE(anyTrack.covariance(eBoundQOverP, eBoundQOverP), 0.0001,
                    1e-6);
  BOOST_CHECK_CLOSE(anyTrack.covariance(eBoundLoc0, eBoundLoc1), 0.0, 1e-6);

  // Verify covariance() map provides reference semantics
  auto covView = anyTrack.covariance();
  covView(eBoundLoc0, eBoundLoc1) = 0.5;
  BOOST_CHECK_CLOSE(track.covariance()(eBoundLoc0, eBoundLoc1), 0.5, 1e-6);

  AnyConstTrackProxy constTrack(track);
  auto constCovView = constTrack.covariance();
  BOOST_CHECK_CLOSE(constCovView(eBoundLoc0, eBoundLoc1), 0.5, 1e-6);
}

BOOST_AUTO_TEST_CASE(AccessParticleHypothesis) {
  VectorTrackContainer vtc;
  VectorMultiTrajectory mtj;
  TrackContainer tc{vtc, mtj};

  auto track = tc.makeTrack();
  fillTestTrack<decltype(tc)>(track);

  AnyMutableTrackProxy anyTrack(track);

  auto ph = anyTrack.particleHypothesis();
  BOOST_CHECK_EQUAL(ph, ParticleHypothesis::pion());
}

BOOST_AUTO_TEST_CASE(AccessDerivedQuantities) {
  VectorTrackContainer vtc;
  VectorMultiTrajectory mtj;
  TrackContainer tc{vtc, mtj};

  auto track = tc.makeTrack();
  fillTestTrack<decltype(tc)>(track);

  AnyMutableTrackProxy anyTrack(track);

  // Check charge
  double expectedCharge =
      ParticleHypothesis::pion().extractCharge(track.qOverP());
  BOOST_CHECK_CLOSE(anyTrack.charge(), expectedCharge, 1e-6);

  // Check momentum
  double expectedMomentum =
      ParticleHypothesis::pion().extractMomentum(track.qOverP());
  BOOST_CHECK_CLOSE(anyTrack.absoluteMomentum(), expectedMomentum, 1e-6);

  // Check transverse momentum
  double expectedPt = std::sin(track.theta()) * expectedMomentum;
  BOOST_CHECK_CLOSE(anyTrack.transverseMomentum(), expectedPt, 1e-6);
}

BOOST_AUTO_TEST_CASE(AccessStatistics) {
  VectorTrackContainer vtc;
  VectorMultiTrajectory mtj;
  TrackContainer tc{vtc, mtj};

  auto track = tc.makeTrack();
  fillTestTrack<decltype(tc)>(track);

  AnyMutableTrackProxy anyTrack(track);

  BOOST_CHECK_EQUAL(anyTrack.chi2(), 12.5f);
  BOOST_CHECK_EQUAL(anyTrack.nDoF(), 10u);
  BOOST_CHECK_EQUAL(anyTrack.nMeasurements(), 8u);
  BOOST_CHECK_EQUAL(anyTrack.nHoles(), 1u);
  BOOST_CHECK_EQUAL(anyTrack.nOutliers(), 0u);
  BOOST_CHECK_EQUAL(anyTrack.nSharedHits(), 2u);
}

BOOST_AUTO_TEST_CASE(AccessTrackStates) {
  VectorTrackContainer vtc;
  VectorMultiTrajectory mtj;
  TrackContainer tc{vtc, mtj};

  auto track = tc.makeTrack();
  fillTestTrack<decltype(tc)>(track);

  // Add some track states
  auto ts1 = tc.trackStateContainer().makeTrackState();
  auto ts2 = tc.trackStateContainer().makeTrackState();
  ts2.previous() = ts1.index();
  auto ts3 = tc.trackStateContainer().makeTrackState();
  ts3.previous() = ts2.index();

  track.tipIndex() = ts3.index();

  AnyMutableTrackProxy anyTrack(track);

  BOOST_CHECK_EQUAL(anyTrack.nTrackStates(), 3u);
}

BOOST_AUTO_TEST_CASE(AccessDynamicColumns) {
  VectorTrackContainer vtc;
  VectorMultiTrajectory mtj;
  TrackContainer tc{vtc, mtj};

  // Add a dynamic column
  tc.addColumn<int>("customInt");
  tc.addColumn<float>("customFloat");

  auto track = tc.makeTrack();
  fillTestTrack<decltype(tc)>(track);

  track.template component<int>("customInt") = 42;
  track.template component<float>("customFloat") = 3.14f;

  AnyMutableTrackProxy anyTrack(track);

  // Check column existence
  BOOST_CHECK(anyTrack.hasColumn("customInt"_hash));
  BOOST_CHECK(anyTrack.hasColumn("customFloat"_hash));
  BOOST_CHECK(!anyTrack.hasColumn("nonExistent"_hash));

  // Access dynamic components using template
  BOOST_CHECK_EQUAL(anyTrack.component<int>("customInt"_hash), 42);
  BOOST_CHECK_CLOSE(anyTrack.component<float>("customFloat"_hash), 3.14f, 1e-6);

  // Also test with compile-time key
  BOOST_CHECK_EQUAL((anyTrack.component<int, "customInt"_hash>()), 42);
  BOOST_CHECK_CLOSE((anyTrack.component<float, "customFloat"_hash>()), 3.14f,
                    1e-6);

  // Mutate via AnyTrack and observe through const handle
  anyTrack.component<int>("customInt"_hash) = 7;
  BOOST_CHECK_EQUAL(anyTrack.component<int>("customInt"_hash), 7);
  AnyConstTrackProxy constTrack(track);
  BOOST_CHECK_EQUAL(constTrack.component<int>("customInt"_hash), 7);
}

BOOST_AUTO_TEST_CASE(ProxyAccessorWithAnyTrack) {
  VectorTrackContainer vtc;
  VectorMultiTrajectory mtj;
  TrackContainer tc{vtc, mtj};

  tc.addColumn<float>("customFloat");
  auto track = tc.makeTrack();
  fillTestTrack<decltype(tc)>(track);
  track.template component<float>("customFloat"_hash) = 1.5f;

  ProxyAccessor<float> mutableAccessor("customFloat");
  ConstProxyAccessor<float> constAccessor("customFloat");

  AnyMutableTrackProxy anyTrack(track);
  BOOST_CHECK_CLOSE(mutableAccessor(anyTrack), 1.5f, 1e-6);
  mutableAccessor(anyTrack) = 2.25f;
  BOOST_CHECK_CLOSE(track.template component<float>("customFloat"_hash), 2.25f,
                    1e-6);

  AnyConstTrackProxy constTrack(track);
  BOOST_CHECK_CLOSE(constAccessor(constTrack), 2.25f, 1e-6);
  BOOST_CHECK(constAccessor.hasColumn(constTrack));
}

BOOST_AUTO_TEST_CASE(TypeErasureHeterogeneousStorage) {
  // Create tracks with different holder types
  VectorTrackContainer vtc1;
  VectorMultiTrajectory mtj1;
  TrackContainer tc1{vtc1, mtj1};

  TrackContainer tc2{VectorTrackContainer{}, VectorMultiTrajectory{}};

  auto vtc3 = std::make_shared<VectorTrackContainer>();
  auto mtj3 = std::make_shared<VectorMultiTrajectory>();
  TrackContainer tc3{vtc3, mtj3};

  auto track1 = tc1.makeTrack();
  auto track2 = tc2.makeTrack();
  auto track3 = tc3.makeTrack();

  fillTestTrack<decltype(tc1)>(track1);
  fillTestTrack<decltype(tc2)>(track2);
  fillTestTrack<decltype(tc3)>(track3);

  track1.chi2() = 10.0f;
  track2.chi2() = 20.0f;
  track3.chi2() = 30.0f;

  // Store in a vector - all have the same type!
  std::vector<AnyMutableTrackProxy> anyTracks;
  anyTracks.emplace_back(track1);
  anyTracks.emplace_back(track2);
  anyTracks.emplace_back(track3);

  BOOST_CHECK_EQUAL(anyTracks.size(), 3u);
  BOOST_CHECK_EQUAL(anyTracks[0].chi2(), 10.0f);
  BOOST_CHECK_EQUAL(anyTracks[1].chi2(), 20.0f);
  BOOST_CHECK_EQUAL(anyTracks[2].chi2(), 30.0f);
}

BOOST_AUTO_TEST_CASE(MemoryFootprint) {
  // Verify that AnyTrack is three pointers (container, handler, index)
  // index is TrackIndexType (std::uint32_t), so on 64-bit it's 2 pointers + 4
  // bytes which gets padded to 3 pointers
  BOOST_CHECK_EQUAL(sizeof(AnyMutableTrackProxy), 3 * sizeof(void*));
  BOOST_CHECK_EQUAL(sizeof(AnyConstTrackProxy), 3 * sizeof(void*));
}

BOOST_AUTO_TEST_CASE(CrossTalkWithTrackProxy) {
  // Test that modifications made through TrackProxy are visible in AnyTrack
  // and vice versa (both reference the same underlying data)
  VectorTrackContainer vtc;
  VectorMultiTrajectory mtj;
  TrackContainer tc{vtc, mtj};

  auto track = tc.makeTrack();
  fillTestTrack<decltype(tc)>(track);

  // Create AnyTrack from TrackProxy
  AnyMutableTrackProxy anyTrack(track);

  // Initial values should match
  BOOST_CHECK_EQUAL(track.chi2(), 12.5f);
  BOOST_CHECK_EQUAL(anyTrack.chi2(), 12.5f);

  // Modify through TrackProxy
  track.chi2() = 42.0f;
  track.nMeasurements() = 15u;
  track.tipIndex() = 99;

  // Changes should be visible in AnyTrack (reads from same container)
  BOOST_CHECK_EQUAL(anyTrack.chi2(), 42.0f);
  BOOST_CHECK_EQUAL(anyTrack.nMeasurements(), 15u);
  BOOST_CHECK_EQUAL(anyTrack.tipIndex(), 99u);

  // Get another TrackProxy to the same track
  auto track2 = tc.getTrack(track.index());

  // Verify both proxies see the same data
  BOOST_CHECK_EQUAL(track2.chi2(), 42.0f);
  BOOST_CHECK_EQUAL(track2.nMeasurements(), 15u);
  BOOST_CHECK_EQUAL(track2.tipIndex(), 99u);

  // Modify through AnyTrack and verify TrackProxy observes the changes
  anyTrack.chi2() = 7.5f;
  anyTrack.nMeasurements() = 5u;
  anyTrack.tipIndex() = 123u;
  anyTrack.parameter(eBoundLoc0) = 5.5;
  anyTrack.covariance(eBoundLoc0, eBoundLoc0) = 0.2;

  BOOST_CHECK_EQUAL(track.chi2(), 7.5f);
  BOOST_CHECK_EQUAL(track.nMeasurements(), 5u);
  BOOST_CHECK_EQUAL(track.tipIndex(), 123u);
  BOOST_CHECK_CLOSE(track.parameters()[eBoundLoc0], 5.5, 1e-6);
  BOOST_CHECK_CLOSE(track.covariance()(eBoundLoc0, eBoundLoc0), 0.2, 1e-6);

  BOOST_CHECK_EQUAL(track2.chi2(), 7.5f);
  BOOST_CHECK_EQUAL(track2.nMeasurements(), 5u);
  BOOST_CHECK_EQUAL(track2.tipIndex(), 123u);
  BOOST_CHECK_CLOSE(track2.parameters()[eBoundLoc0], 5.5, 1e-6);
  BOOST_CHECK_CLOSE(track2.covariance()(eBoundLoc0, eBoundLoc0), 0.2, 1e-6);

  // Modify through the second proxy
  track2.chi2() = 100.0f;

  // Changes should be visible everywhere
  BOOST_CHECK_EQUAL(track.chi2(), 100.0f);
  BOOST_CHECK_EQUAL(anyTrack.chi2(), 100.0f);
}

BOOST_AUTO_TEST_CASE(ConstCorrectnessCrossTalk) {
  // Test that const AnyTrack can be created from mutable TrackProxy
  // and both see the same data
  VectorTrackContainer vtc;
  VectorMultiTrajectory mtj;
  TrackContainer tc{vtc, mtj};

  auto track = tc.makeTrack();
  fillTestTrack<decltype(tc)>(track);
  track.chi2() = 25.0f;

  // Create const AnyTrack from mutable TrackProxy
  AnyConstTrackProxy constAnyTrack(track);

  // Both should see the same data
  BOOST_CHECK_EQUAL(track.chi2(), 25.0f);
  BOOST_CHECK_EQUAL(constAnyTrack.chi2(), 25.0f);

  // Modify through mutable proxy
  track.chi2() = 50.0f;

  // Const AnyTrack should see the change
  BOOST_CHECK_EQUAL(constAnyTrack.chi2(), 50.0f);
}

BOOST_AUTO_TEST_CASE(ConstCorrectnessRestrictions) {
  // Test const-correctness restrictions on constructor
  VectorTrackContainer vtc;
  VectorMultiTrajectory mtj;
  TrackContainer tc{vtc, mtj};

  auto mutableTrack = tc.makeTrack();
  fillTestTrack<decltype(tc)>(mutableTrack);

  auto constTrack = tc.getTrack(mutableTrack.index());

  // These should compile:
  // 1. AnyMutableTrackProxy from mutable track proxy
  AnyMutableTrackProxy anyMutable1(mutableTrack);

  // 2. AnyConstTrackProxy from mutable track proxy
  AnyConstTrackProxy anyConst1(mutableTrack);

  // 3. AnyConstTrackProxy from const track proxy
  AnyConstTrackProxy anyConst2(constTrack);

  // The following should NOT compile (uncomment to test):
  // AnyMutableTrackProxy anyMutable2(constTrack);  // Error: cannot create
  // mutable from const

  // Verify all see the same data
  mutableTrack.chi2() = 77.0f;
  BOOST_CHECK_EQUAL(anyMutable1.chi2(), 77.0f);
  BOOST_CHECK_EQUAL(anyConst1.chi2(), 77.0f);
  BOOST_CHECK_EQUAL(anyConst2.chi2(), 77.0f);
}

BOOST_AUTO_TEST_SUITE_END()

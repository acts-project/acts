// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Definitions/Units.hpp"
#include "Acts/EventData/ProxyAccessor.hpp"
#include "Acts/EventData/VectorMultiTrajectory.hpp"
#include "Acts/EventData/VectorTrackContainer.hpp"
#include "Acts/Surfaces/CurvilinearSurface.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/Utilities/Diagnostics.hpp"
#include "Acts/Utilities/Zip.hpp"

#include <algorithm>
#include <numeric>

using namespace Acts;
using namespace Acts::HashedStringLiteral;
using IndexType = TrackIndexType;
constexpr auto kInvalid = kTrackIndexInvalid;
using namespace Acts::UnitLiterals;

template <typename track_container_t, typename traj_t,
          template <typename> class holder_t>
struct Factory {};

template <typename track_container_t, typename traj_t>
struct Factory<track_container_t, traj_t, RefHolder> {
  using track_container_type =
      TrackContainer<track_container_t, traj_t, RefHolder>;

  track_container_t vtc;
  traj_t mtj;
  track_container_type tc{vtc, mtj};

  auto& trackContainer() { return tc; }
  auto& trackStateContainer() { return mtj; }
  auto& backend() { return vtc; }
};

template <typename track_container_t, typename traj_t>
struct Factory<track_container_t, traj_t, ValueHolder> {
  using track_container_type =
      TrackContainer<track_container_t, traj_t, ValueHolder>;

  track_container_type tc{track_container_t{}, traj_t{}};

  auto& trackContainer() { return tc; }
  auto& trackStateContainer() { return tc.trackStateContainer(); }
  auto& backend() { return tc.container(); }
};

template <typename track_container_t, typename traj_t>
struct Factory<track_container_t, traj_t, std::shared_ptr> {
  using track_container_type =
      TrackContainer<track_container_t, traj_t, std::shared_ptr>;

  std::shared_ptr<track_container_t> vtc{std::make_shared<track_container_t>()};
  std::shared_ptr<traj_t> mtj{std::make_shared<traj_t>()};
  track_container_type tc{vtc, mtj};

  auto& trackContainer() { return tc; }
  auto& trackStateContainer() { return *mtj; }
  auto& backend() { return *vtc; }
};

template <typename track_container_t, typename traj_t,
          template <typename> class... holders>
using holder_types_t =
    std::tuple<Factory<track_container_t, traj_t, holders>...>;

using holder_types = holder_types_t<VectorTrackContainer, VectorMultiTrajectory,
                                    // detail_tc::ValueHolder,
                                    // detail_tc::RefHolder,
                                    std::shared_ptr>;

using const_holder_types =
    holder_types_t<ConstVectorTrackContainer, ConstVectorMultiTrajectory,
                   ValueHolder, RefHolder, std::shared_ptr>;

namespace ActsTests {

BOOST_AUTO_TEST_SUITE(EventDataSuite)

BOOST_AUTO_TEST_CASE(BuildDefaultHolder) {
  VectorMultiTrajectory mtj{};
  VectorTrackContainer vtc{};
  TrackContainer tc{vtc, mtj};

  static_assert(
      std::is_same_v<decltype(tc),
                     TrackContainer<VectorTrackContainer, VectorMultiTrajectory,
                                    RefHolder>>,
      "Incorrect deduced type");
  BOOST_CHECK_EQUAL(&mtj, &tc.trackStateContainer());
  BOOST_CHECK_EQUAL(&vtc, &tc.container());
  tc.addTrack();

  std::decay_t<decltype(tc)> copy = tc;
  BOOST_CHECK_EQUAL(&mtj, &copy.trackStateContainer());
  BOOST_CHECK_EQUAL(&vtc, &copy.container());
}

BOOST_AUTO_TEST_CASE(BuildValueHolder) {
  {
    VectorMultiTrajectory mtj{};
    VectorTrackContainer vtc{};
    TrackContainer tc{std::move(vtc), std::move(mtj)};
    static_assert(
        std::is_same_v<decltype(tc),
                       TrackContainer<VectorTrackContainer,
                                      VectorMultiTrajectory, ValueHolder>>,
        "Incorrect deduced type");
    std::decay_t<decltype(tc)> copy = tc;
    BOOST_CHECK_NE(&tc.trackStateContainer(), &copy.trackStateContainer());
    BOOST_CHECK_NE(&tc.container(), &copy.container());
  }
  {
    TrackContainer tc{VectorTrackContainer{}, VectorMultiTrajectory{}};

    static_assert(
        std::is_same_v<decltype(tc),
                       TrackContainer<VectorTrackContainer,
                                      VectorMultiTrajectory, ValueHolder>>,
        "Incorrect deduced type");
    tc.addTrack();
    std::decay_t<decltype(tc)> copy = tc;
    BOOST_CHECK_NE(&tc.trackStateContainer(), &copy.trackStateContainer());
    BOOST_CHECK_NE(&tc.container(), &copy.container());
  }
}

BOOST_AUTO_TEST_CASE(BuildRefHolder) {
  VectorMultiTrajectory mtj{};
  VectorTrackContainer vtc{};
  TrackContainer<VectorTrackContainer, VectorMultiTrajectory, RefHolder> tc{
      vtc, mtj};

  static_assert(
      std::is_same_v<decltype(tc),
                     TrackContainer<VectorTrackContainer, VectorMultiTrajectory,
                                    RefHolder>>,
      "Incorrect deduced type");
  BOOST_CHECK_EQUAL(&mtj, &tc.trackStateContainer());
  BOOST_CHECK_EQUAL(&vtc, &tc.container());
  tc.addTrack();
  std::decay_t<decltype(tc)> copy = tc;
  BOOST_CHECK_EQUAL(&mtj, &copy.trackStateContainer());
  BOOST_CHECK_EQUAL(&vtc, &copy.container());
}

BOOST_AUTO_TEST_CASE(BuildSharedPtr) {
  auto mtj = std::make_shared<VectorMultiTrajectory>();
  auto vtc = std::make_shared<VectorTrackContainer>();
  TrackContainer<VectorTrackContainer, VectorMultiTrajectory, std::shared_ptr>
      tc{vtc, mtj};

  static_assert(
      std::is_same_v<decltype(tc),
                     TrackContainer<VectorTrackContainer, VectorMultiTrajectory,
                                    std::shared_ptr>>,
      "Incorrect deduced type");
  BOOST_CHECK_EQUAL(mtj.get(), &tc.trackStateContainer());
  BOOST_CHECK_EQUAL(vtc.get(), &tc.container());
  tc.addTrack();
  std::decay_t<decltype(tc)> copy = tc;
  BOOST_CHECK_EQUAL(mtj.get(), &copy.trackStateContainer());
  BOOST_CHECK_EQUAL(vtc.get(), &copy.container());
}

BOOST_AUTO_TEST_CASE(BuildConvenience) {
  VectorMultiTrajectory mtj{};
  VectorTrackContainer vtc{};
  TrackContainer tc{vtc, mtj};

  BOOST_CHECK_EQUAL(tc.size(), 0);
  auto track1 = tc.makeTrack();
  BOOST_CHECK_EQUAL(tc.size(), 1);
  auto track2 = tc.makeTrack();
  BOOST_CHECK_EQUAL(tc.size(), 2);

  BOOST_CHECK_EQUAL(track1.index(), 0);
  BOOST_CHECK_EQUAL(track2.index(), 1);
}

BOOST_AUTO_TEST_CASE_TEMPLATE(Build, factory_t, holder_types) {
  factory_t factory;

  auto& tc = factory.trackContainer();

  static_assert(std::is_same_v<std::decay_t<decltype(tc)>,
                               typename factory_t::track_container_type>,
                "Incorrect deduction");

  static_assert(!std::decay_t<decltype(tc)>::ReadOnly,
                "Should not be read only");
  BOOST_CHECK(!tc.ReadOnly);

  auto idx = tc.addTrack();
  auto t = tc.getTrack(idx);
  auto t2 = tc.getTrack(idx);
  t.template component<IndexType, "tipIndex"_hash>() = 5;

  BOOST_CHECK_EQUAL((t.template component<IndexType, "tipIndex"_hash>()), 5);
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

  std::shared_ptr<PlaneSurface> surface =
      CurvilinearSurface(Vector3{-3_m, 0., 0.}, Vector3{1., 0., 0})
          .planeSurface();

  t.setReferenceSurface(surface);
  BOOST_CHECK_EQUAL(surface.get(), &t.referenceSurface());

  ProxyAccessor<unsigned int> accNMeasuements("nMeasurements");
  ConstProxyAccessor<unsigned int> caccNMeasuements("nMeasurements");

  t.nMeasurements() = 42;
  BOOST_CHECK_EQUAL(t2.nMeasurements(), 42);
  BOOST_CHECK_EQUAL(accNMeasuements(t), 42);
  accNMeasuements(t) = 89;
  BOOST_CHECK_EQUAL(t2.nMeasurements(), 89);
  BOOST_CHECK_EQUAL(caccNMeasuements(t), 89);

  // Test hasColumn functionality
  BOOST_CHECK(accNMeasuements.hasColumn(t));
  BOOST_CHECK(caccNMeasuements.hasColumn(t));

  ProxyAccessor<std::string> accNonExistent("nonExistentColumn");
  ConstProxyAccessor<std::string> caccNonExistent("nonExistentColumn");
  BOOST_CHECK(!accNonExistent.hasColumn(t));
  BOOST_CHECK(!caccNonExistent.hasColumn(t));

  // does not compile
  // caccNMeasuements(t) = 66;

  t2.nHoles() = 67;
  BOOST_CHECK_EQUAL(t.nHoles(), 67);

  t2.nOutliers() = 68;
  BOOST_CHECK_EQUAL(t.nOutliers(), 68);

  t2.nSharedHits() = 69;
  BOOST_CHECK_EQUAL(t.nSharedHits(), 69);

  t2.chi2() = 555.0;
  BOOST_CHECK_EQUAL(t2.chi2(), 555.0);

  t2.nDoF() = 123;
  BOOST_CHECK_EQUAL(t2.nDoF(), 123);

  // const checks: should not compile
  // const auto& ctc = tc;
  // ctc.getTrack(idx).covariance().setRandom();
  // const auto& ctp = t;
  // ctp.covariance().setRandom();
}

BOOST_AUTO_TEST_CASE(CopyTracksIncludingDynamicColumns) {
  // mutable source
  VectorTrackContainer vtc{};
  VectorMultiTrajectory mtj{};
  TrackContainer tc{vtc, mtj};
  tc.addColumn<std::size_t>("counter");
  tc.addColumn<bool>("odd");
  mtj.addColumn<std::size_t>("ts_counter");
  mtj.addColumn<bool>("ts_odd");

  TrackContainer tc2{VectorTrackContainer{}, VectorMultiTrajectory{}};
  // doesn't have the dynamic column

  VectorTrackContainer vtc3{};
  VectorMultiTrajectory mtj3{};
  mtj3.addColumn<std::size_t>("ts_counter");
  mtj3.addColumn<bool>("ts_odd");

  TrackContainer tc3{vtc3, mtj3};

  tc3.addColumn<std::size_t>("counter");
  tc3.addColumn<bool>("odd");

  for (std::size_t i = 0; i < 10; i++) {
    auto t = tc.makeTrack();
    auto ts = t.appendTrackState();
    ts.predicted() = BoundVector::Ones();
    ts.component<std::size_t, "ts_counter"_hash>() = i;

    ts = t.appendTrackState();
    ts.predicted().setOnes();
    ts.predicted() *= 2;
    ts.component<std::size_t, "ts_counter"_hash>() = i + 1;

    ts = t.appendTrackState();
    ts.predicted().setOnes();
    ts.predicted() *= 3;
    ts.component<std::size_t, "ts_counter"_hash>() = i + 2;

    t.template component<std::size_t>("counter") = i;
    t.template component<bool>("odd") = i % 2 == 0;

    auto t2 = tc2.makeTrack();
    BOOST_CHECK_THROW(t2.copyFrom(t),
                      std::invalid_argument);  // this should fail

    auto t3 = tc3.makeTrack();
    t3.copyFrom(t);  // this should work

    BOOST_CHECK_NE(t3.tipIndex(), kInvalid);
    BOOST_CHECK_GT(t3.nTrackStates(), 0);
    BOOST_REQUIRE_EQUAL(t.nTrackStates(), t3.nTrackStates());

    for (auto [tsa, tsb] :
         zip(t.trackStatesReversed(), t3.trackStatesReversed())) {
      BOOST_CHECK_EQUAL(tsa.predicted(), tsb.predicted());

      BOOST_CHECK_EQUAL(
          (tsa.template component<std::size_t, "ts_counter"_hash>()),
          (tsb.template component<std::size_t, "ts_counter"_hash>()));

      BOOST_CHECK_EQUAL((tsa.template component<bool, "ts_odd"_hash>()),
                        (tsb.template component<bool, "ts_odd"_hash>()));
    }

    BOOST_CHECK_EQUAL(t.template component<std::size_t>("counter"),
                      t3.template component<std::size_t>("counter"));
    BOOST_CHECK_EQUAL(t.template component<bool>("odd"),
                      t3.template component<bool>("odd"));
  }

  std::size_t before = mtj.size();
  TrackContainer tc4{ConstVectorTrackContainer{vtc},
                     ConstVectorMultiTrajectory{mtj}};

  BOOST_REQUIRE_EQUAL(tc4.trackStateContainer().size(), before);

  VectorTrackContainer vtc5{};
  VectorMultiTrajectory mtj5{};
  mtj5.addColumn<std::size_t>("ts_counter");
  mtj5.addColumn<bool>("ts_odd");

  TrackContainer tc5{vtc5, mtj5};
  tc5.addColumn<std::size_t>("counter");
  tc5.addColumn<bool>("odd");

  for (std::size_t i = 0; i < 10; i++) {
    auto t4 = tc4.getTrack(i);  // const source!
    BOOST_CHECK_NE(t4.nTrackStates(), 0);

    auto t5 = tc5.makeTrack();
    t5.copyFrom(t4);  // this should work

    BOOST_CHECK_NE(t5.tipIndex(), kInvalid);
    BOOST_CHECK_GT(t5.nTrackStates(), 0);
    BOOST_REQUIRE_EQUAL(t4.nTrackStates(), t5.nTrackStates());

    for (auto [tsa, tsb] :
         zip(t4.trackStatesReversed(), t5.trackStatesReversed())) {
      BOOST_CHECK_EQUAL(tsa.predicted(), tsb.predicted());
    }

    BOOST_CHECK_EQUAL(t4.template component<std::size_t>("counter"),
                      t5.template component<std::size_t>("counter"));
    BOOST_CHECK_EQUAL(t4.template component<bool>("odd"),
                      t5.template component<bool>("odd"));
  }
}

BOOST_AUTO_TEST_CASE(ReverseTrackStates) {
  VectorTrackContainer vtc{};
  VectorMultiTrajectory mtj{};
  TrackContainer tc{vtc, mtj};

  auto t = tc.makeTrack();

  for (std::size_t i = 0; i < 4; i++) {
    auto ts = t.appendTrackState();
    ts.jacobian() = BoundMatrix::Identity() * i;
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
    BOOST_CHECK_EQUAL(ts.jacobian(), BoundMatrix::Identity() * e);
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
    BOOST_CHECK_EQUAL(ts.jacobian(), BoundMatrix::Identity() * e);
  }

  // back to original!
  t.reverseTrackStates();

  // jacobians stay with their track states
  for (const auto [e, ts] : zip(exp, t.trackStates())) {
    BOOST_CHECK_EQUAL(ts.jacobian(), BoundMatrix::Identity() * e);
  }

  // reverse with jacobians
  t.reverseTrackStates(true);

  std::ranges::reverse(exp);
  std::rotate(exp.rbegin(), std::next(exp.rbegin()), exp.rend());

  for (const auto [e, ts] : zip(exp, t.trackStates())) {
    Acts::BoundMatrix expJac;
    if (e == 0) {
      expJac = BoundMatrix::Zero();
    } else {
      expJac = (BoundMatrix::Identity() * e).inverse();
    }

    BOOST_CHECK_EQUAL(ts.jacobian(), expJac);
  }

  // now back to original order, revert jacobians again
  t.reverseTrackStates(true);

  // reset exp to range(0, N)
  std::iota(exp.begin(), exp.end(), 0);

  for (const auto [e, ts] : zip(exp, t.trackStates())) {
    BOOST_CHECK_EQUAL(ts.jacobian(), BoundMatrix::Identity() * e);
  }
}

BOOST_AUTO_TEST_CASE(CopyTrackProxyCalibrated) {
  VectorTrackContainer vtc{};
  VectorMultiTrajectory mtj{};
  TrackContainer tc{vtc, mtj};

  constexpr static std::size_t kMeasurementSize = 3;

  auto track1 = tc.makeTrack();
  auto ts = track1.appendTrackState(TrackStatePropMask::Calibrated);
  ts.allocateCalibrated(kMeasurementSize);
  ts.setProjectorSubspaceIndices(BoundSubspaceIndices{});

  auto tsCopy = track1.appendTrackState(TrackStatePropMask::Calibrated);
  tsCopy.copyFrom(ts, TrackStatePropMask::Calibrated, false);

  BOOST_CHECK_EQUAL(ts.calibratedSize(), tsCopy.calibratedSize());
}

BOOST_AUTO_TEST_CASE(ProxyAccessorHasColumn) {
  VectorTrackContainer vtc{};
  VectorMultiTrajectory mtj{};
  TrackContainer tc{vtc, mtj};

  // Add a custom column to the track state container
  mtj.addColumn<float>("customFloat");
  // Add a custom column to the track container
  tc.addColumn<std::string>("customString");

  auto track = tc.makeTrack();
  auto ts = track.appendTrackState(TrackStatePropMask::Predicted);
  ts.predicted() = BoundVector::Zero();

  // Set values in the custom columns
  ts.component<float>("customFloat") = 42.5f;
  track.component<std::string>("customString") = "test";

  // Test ALL known track container columns
  // From VectorTrackContainer.hpp hasColumn_impl
  BOOST_CHECK(ConstProxyAccessor<IndexType>("tipIndex").hasColumn(track));
  BOOST_CHECK(ConstProxyAccessor<IndexType>("stemIndex").hasColumn(track));
  BOOST_CHECK(ConstProxyAccessor<BoundVector>("params").hasColumn(track));
  BOOST_CHECK(ConstProxyAccessor<BoundMatrix>("cov").hasColumn(track));
  BOOST_CHECK(
      ConstProxyAccessor<unsigned int>("nMeasurements").hasColumn(track));
  BOOST_CHECK(ConstProxyAccessor<unsigned int>("nHoles").hasColumn(track));
  BOOST_CHECK(ConstProxyAccessor<float>("chi2").hasColumn(track));
  BOOST_CHECK(ConstProxyAccessor<unsigned int>("ndf").hasColumn(track));
  BOOST_CHECK(ConstProxyAccessor<unsigned int>("nOutliers").hasColumn(track));
  BOOST_CHECK(ConstProxyAccessor<unsigned int>("nSharedHits").hasColumn(track));

  // Test custom track container column
  BOOST_CHECK(ConstProxyAccessor<std::string>("customString").hasColumn(track));

  // Test ALL known track state columns
  // From VectorMultiTrajectory.hpp hasColumn_impl
  BOOST_CHECK(ConstProxyAccessor<IndexType>("previous").hasColumn(ts));
  BOOST_CHECK(ConstProxyAccessor<IndexType>("next").hasColumn(ts));
  BOOST_CHECK(ConstProxyAccessor<IndexType>("predicted").hasColumn(ts));
  BOOST_CHECK(ConstProxyAccessor<IndexType>("filtered").hasColumn(ts));
  BOOST_CHECK(ConstProxyAccessor<IndexType>("smoothed").hasColumn(ts));
  BOOST_CHECK(ConstProxyAccessor<IndexType>("calibrated").hasColumn(ts));
  BOOST_CHECK(ConstProxyAccessor<IndexType>("calibratedCov").hasColumn(ts));
  BOOST_CHECK(ConstProxyAccessor<IndexType>("jacobian").hasColumn(ts));
  BOOST_CHECK(ConstProxyAccessor<IndexType>("projector").hasColumn(ts));
  BOOST_CHECK(
      ConstProxyAccessor<std::optional<SourceLink>>("uncalibratedSourceLink")
          .hasColumn(ts));
  BOOST_CHECK(
      ConstProxyAccessor<std::shared_ptr<const Surface>>("referenceSurface")
          .hasColumn(ts));
  BOOST_CHECK(ConstProxyAccessor<std::uint8_t>("measdim").hasColumn(ts));
  BOOST_CHECK(ConstProxyAccessor<float>("chi2").hasColumn(ts));
  BOOST_CHECK(ConstProxyAccessor<double>("pathLength").hasColumn(ts));
  BOOST_CHECK(
      ConstProxyAccessor<TrackStateType::raw_type>("typeFlags").hasColumn(ts));

  // Test custom track state column
  BOOST_CHECK(ConstProxyAccessor<float>("customFloat").hasColumn(ts));

  // Test with non-existent columns
  BOOST_CHECK(!ConstProxyAccessor<int>("nonExistentColumn").hasColumn(ts));
  BOOST_CHECK(!ConstProxyAccessor<std::string>("nonExistentTrackColumn")
                   .hasColumn(track));

  // Test that we can actually access the custom columns
  ProxyAccessor<float> accCustomFloat("customFloat");
  ConstProxyAccessor<float> caccCustomFloat("customFloat");
  BOOST_CHECK_EQUAL(accCustomFloat(ts), 42.5f);
  BOOST_CHECK_EQUAL(caccCustomFloat(ts), 42.5f);

  ProxyAccessor<std::string> accCustomString("customString");
  ConstProxyAccessor<std::string> caccCustomString("customString");
  BOOST_CHECK_EQUAL(accCustomString(track), "test");
  BOOST_CHECK_EQUAL(caccCustomString(track), "test");
}

consteval ProxyAccessor<int> makeProxyAccessor() {
  return ProxyAccessor<int>("test_consteval");
}

BOOST_AUTO_TEST_CASE(ProxyAccessorConstexprConstruction) {
  static_assert(ProxyAccessor<int>("test").key == "test"_hash,
                "ProxyAccessor should be constructible with a string key");

  constexpr auto acc = makeProxyAccessor();

  static_assert(acc.key == "test_consteval"_hash,
                "ProxyAccessor should be constructible at compile time");
}

BOOST_AUTO_TEST_CASE(CopyFromWithoutStatesInvalidatesIndices) {
  VectorTrackContainer vtc{};
  VectorMultiTrajectory mtj{};
  TrackContainer tc{vtc, mtj};

  // Create source track with track states
  auto sourceTrack = tc.makeTrack();
  sourceTrack.appendTrackState();
  sourceTrack.appendTrackState();
  sourceTrack.appendTrackState();
  sourceTrack.linkForward();

  // Verify source track has valid indices
  BOOST_CHECK_NE(sourceTrack.tipIndex(), kInvalid);
  BOOST_CHECK_NE(sourceTrack.stemIndex(), kInvalid);
  BOOST_CHECK_EQUAL(sourceTrack.nTrackStates(), 3);

  // Set some track properties
  sourceTrack.nMeasurements() = 42;
  sourceTrack.nHoles() = 5;
  sourceTrack.chi2() = 123.45f;

  // Create destination track that also has track states initially
  auto destTrack = tc.makeTrack();
  destTrack.appendTrackState();
  destTrack.appendTrackState();
  destTrack.linkForward();

  // Verify destination has valid indices before copy
  BOOST_CHECK_NE(destTrack.tipIndex(), kInvalid);
  BOOST_CHECK_NE(destTrack.stemIndex(), kInvalid);
  BOOST_CHECK_EQUAL(destTrack.nTrackStates(), 2);

  // Copy without states
  destTrack.copyFromWithoutStates(sourceTrack);

  // Verify that tip and stem indices are now invalid
  BOOST_CHECK_EQUAL(destTrack.tipIndex(), kInvalid);
  BOOST_CHECK_EQUAL(destTrack.stemIndex(), kInvalid);

  // Verify that track properties were copied
  BOOST_CHECK_EQUAL(destTrack.nMeasurements(), 42);
  BOOST_CHECK_EQUAL(destTrack.nHoles(), 5);
  BOOST_CHECK_EQUAL(destTrack.chi2(), 123.45f);

  // Verify that nTrackStates returns 0 since indices are invalid
  BOOST_CHECK_EQUAL(destTrack.nTrackStates(), 0);

  // Verify that the original track states are still in the container
  // but not accessible through the destination track
  BOOST_CHECK_EQUAL(sourceTrack.nTrackStates(), 3);
  BOOST_CHECK_NE(sourceTrack.tipIndex(), kInvalid);
}

BOOST_AUTO_TEST_CASE(CopyFromDeepCopyFunctionality) {
  VectorTrackContainer vtc{};
  VectorMultiTrajectory mtj{};
  TrackContainer tc{vtc, mtj};

  // Create source track with track states and set unique data
  auto sourceTrack = tc.makeTrack();

  // Set track-level properties
  sourceTrack.nMeasurements() = 15;
  sourceTrack.nHoles() = 3;
  sourceTrack.nOutliers() = 2;
  sourceTrack.chi2() = 42.5f;
  sourceTrack.nDoF() = 10;

  // Create track states with distinctive predicted parameters
  auto ts1 = sourceTrack.appendTrackState();
  ts1.predicted() = BoundVector::Ones() * 1.0;

  auto ts2 = sourceTrack.appendTrackState();
  ts2.predicted() = BoundVector::Ones() * 2.0;

  auto ts3 = sourceTrack.appendTrackState();
  ts3.predicted() = BoundVector::Ones() * 3.0;

  // Verify source track setup
  BOOST_CHECK_EQUAL(sourceTrack.nTrackStates(), 3);

  // Create destination track and perform deep copy
  auto destTrack = tc.makeTrack();
  destTrack.copyFrom(sourceTrack);

  // Verify track-level properties were copied
  BOOST_CHECK_EQUAL(destTrack.nMeasurements(), 15);
  BOOST_CHECK_EQUAL(destTrack.nHoles(), 3);
  BOOST_CHECK_EQUAL(destTrack.nOutliers(), 2);
  BOOST_CHECK_EQUAL(destTrack.chi2(), 42.5f);
  BOOST_CHECK_EQUAL(destTrack.nDoF(), 10);

  // Verify track state structure
  BOOST_CHECK_EQUAL(destTrack.nTrackStates(), 3);
  BOOST_CHECK_NE(destTrack.tipIndex(), kInvalid);
  BOOST_CHECK_NE(destTrack.stemIndex(), kInvalid);

  // Verify track is forward-linked (can iterate forward)
  BOOST_CHECK(destTrack.innermostTrackState().has_value());

  // Verify track state data was copied correctly (order preserved)
  std::vector<BoundVector> sourceParams;
  for (const auto& ts : sourceTrack.trackStatesReversed()) {
    sourceParams.insert(sourceParams.begin(), ts.predicted());
  }

  std::vector<BoundVector> destParams;
  for (const auto& ts : destTrack.trackStatesReversed()) {
    destParams.insert(destParams.begin(), ts.predicted());
  }

  BOOST_REQUIRE_EQUAL(sourceParams.size(), destParams.size());
  for (std::size_t i = 0; i < sourceParams.size(); ++i) {
    BOOST_CHECK_EQUAL(sourceParams[i], destParams[i]);
  }

  // Verify track states have different indices (new track states were created)
  std::vector<IndexType> sourceIndices;
  for (const auto& ts : sourceTrack.trackStatesReversed()) {
    sourceIndices.push_back(ts.index());
  }

  std::vector<IndexType> destIndices;
  for (const auto& ts : destTrack.trackStatesReversed()) {
    destIndices.push_back(ts.index());
  }

  BOOST_REQUIRE_EQUAL(sourceIndices.size(), destIndices.size());
  for (std::size_t i = 0; i < sourceIndices.size(); ++i) {
    BOOST_CHECK_NE(sourceIndices[i], destIndices[i]);
  }

  // Verify forward iteration works (result of forward-linking)
  std::vector<BoundVector> forwardParams;
  for (const auto& ts : destTrack.trackStates()) {
    forwardParams.push_back(ts.predicted());
  }

  // Forward iteration should give same parameters in same order
  BOOST_REQUIRE_EQUAL(destParams.size(), forwardParams.size());
  for (std::size_t i = 0; i < destParams.size(); ++i) {
    BOOST_CHECK_EQUAL(destParams[i], forwardParams[i]);
  }

  // Verify tracks are independent (modifying dest doesn't affect source)
  destTrack.nMeasurements() = 99;
  BOOST_CHECK_EQUAL(sourceTrack.nMeasurements(), 15);
  BOOST_CHECK_EQUAL(destTrack.nMeasurements(), 99);
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace ActsTests

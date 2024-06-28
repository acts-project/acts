// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/Definitions/Units.hpp"
#include "Acts/EventData/MultiTrajectory.hpp"
#include "Acts/EventData/ProxyAccessor.hpp"
#include "Acts/EventData/TrackContainer.hpp"
#include "Acts/EventData/TrackHelpers.hpp"
#include "Acts/EventData/TrackProxy.hpp"
#include "Acts/EventData/TrackStatePropMask.hpp"
#include "Acts/EventData/VectorMultiTrajectory.hpp"
#include "Acts/EventData/VectorTrackContainer.hpp"
#include "Acts/EventData/detail/GenerateParameters.hpp"
#include "Acts/EventData/detail/TestTrackState.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/HashedString.hpp"
#include "Acts/Utilities/Holders.hpp"
#include "Acts/Utilities/Zip.hpp"

#include <cstddef>
#include <iterator>
#include <memory>
#include <numeric>
#include <random>
#include <stdexcept>
#include <tuple>
#include <type_traits>
#include <utility>
#include <vector>

namespace {

using namespace Acts::UnitLiterals;

using namespace Acts;
using namespace Acts::HashedStringLiteral;
using namespace Acts::detail::Test;

using MultiTrajectoryTraits::IndexType;

const GeometryContext gctx;
// fixed seed for reproducible tests
std::default_random_engine rng(31415);

// template <template <typename> class holder_t>
// using track_container_t =
// TrackContainer<VectorTrackContainer, VectorMultiTrajectory, holder_t>;

template <typename track_container_t, typename traj_t,
          template <typename> class holder_t>
struct Factory {};

template <typename track_container_t, typename traj_t>
struct Factory<track_container_t, traj_t, detail::RefHolder> {
  using track_container_type =
      TrackContainer<track_container_t, traj_t, detail::RefHolder>;

  track_container_t vtc;
  traj_t mtj;
  track_container_type tc{vtc, mtj};

  auto& trackContainer() { return tc; }
  auto& trackStateContainer() { return mtj; }
  auto& backend() { return vtc; }
};

template <typename track_container_t, typename traj_t>
struct Factory<track_container_t, traj_t, detail::ValueHolder> {
  using track_container_type =
      TrackContainer<track_container_t, traj_t, detail::ValueHolder>;

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
                   detail::ValueHolder, detail::RefHolder, std::shared_ptr>;

}  // namespace

BOOST_AUTO_TEST_SUITE(EventDataTrack)

BOOST_AUTO_TEST_CASE(BuildDefaultHolder) {
  VectorMultiTrajectory mtj{};
  VectorTrackContainer vtc{};
  TrackContainer tc{vtc, mtj};

  static_assert(
      std::is_same_v<decltype(tc),
                     TrackContainer<VectorTrackContainer, VectorMultiTrajectory,
                                    detail::RefHolder>>,
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
        std::is_same_v<decltype(tc), TrackContainer<VectorTrackContainer,
                                                    VectorMultiTrajectory,
                                                    detail::ValueHolder>>,
        "Incorrect deduced type");
    std::decay_t<decltype(tc)> copy = tc;
    BOOST_CHECK_NE(&tc.trackStateContainer(), &copy.trackStateContainer());
    BOOST_CHECK_NE(&tc.container(), &copy.container());
  }
  {
    TrackContainer tc{VectorTrackContainer{}, VectorMultiTrajectory{}};

    static_assert(
        std::is_same_v<decltype(tc), TrackContainer<VectorTrackContainer,
                                                    VectorMultiTrajectory,
                                                    detail::ValueHolder>>,
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
  TrackContainer<VectorTrackContainer, VectorMultiTrajectory, detail::RefHolder>
      tc{vtc, mtj};

  static_assert(
      std::is_same_v<decltype(tc),
                     TrackContainer<VectorTrackContainer, VectorMultiTrajectory,
                                    detail::RefHolder>>,
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

  auto surface = Acts::Surface::makeShared<Acts::PlaneSurface>(
      Acts::Vector3{-3_m, 0., 0.}, Acts::Vector3{1., 0., 0});

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

BOOST_AUTO_TEST_CASE_TEMPLATE(TrackStateAccess, factory_t, holder_types) {
  factory_t factory;
  auto& tc = factory.trackContainer();

  VectorMultiTrajectory& traj = factory.trackStateContainer();

  auto mkts = [&](auto prev) {
    if constexpr (std::is_same_v<decltype(prev), IndexType>) {
      auto ts = traj.makeTrackState(TrackStatePropMask::All, prev);
      TestTrackState pc(rng, 2u);
      fillTrackState<VectorMultiTrajectory>(pc, TrackStatePropMask::All, ts);
      return ts;
    } else {
      auto ts = traj.makeTrackState(TrackStatePropMask::All, prev.index());
      TestTrackState pc(rng, 2u);
      fillTrackState<VectorMultiTrajectory>(pc, TrackStatePropMask::All, ts);
      return ts;
    }
  };

  auto ts1 = mkts(MultiTrajectoryTraits::kInvalid);
  auto ts2 = mkts(ts1);
  auto ts3 = mkts(ts2);
  auto ts4 = mkts(ts3);
  auto ts5 = mkts(ts4);

  auto t = tc.makeTrack();
  t.tipIndex() = ts5.index();

  std::vector<IndexType> act;
  for (const auto& ts : t.trackStatesReversed()) {
    act.push_back(ts.index());
  }

  std::vector<IndexType> exp;
  exp.resize(5);
  std::iota(exp.rbegin(), exp.rend(), 0);
  BOOST_CHECK_EQUAL_COLLECTIONS(act.begin(), act.end(), exp.begin(), exp.end());

  const auto& ct = t;

  for (const auto& ts : ct.trackStatesReversed()) {
    (void)ts;
  }

  BOOST_CHECK_EQUAL(t.nTrackStates(), 5);

  auto tNone = tc.makeTrack();
  BOOST_CHECK_EQUAL(tNone.nTrackStates(), 0);

  auto tsRange = tNone.trackStatesReversed();
  BOOST_CHECK(tsRange.begin() == tsRange.end());

  std::size_t i = 0;
  for (const auto& state : tNone.trackStatesReversed()) {
    (void)state;
    i++;
  }
  BOOST_CHECK_EQUAL(i, 0);
}

BOOST_AUTO_TEST_CASE_TEMPLATE(TrackIterator, factory_t, holder_types) {
  factory_t factory;
  auto& tc = factory.trackContainer();

  for (unsigned int i = 0; i < 10; i++) {
    auto t = tc.makeTrack();
    t.tipIndex() = i;
  }
  BOOST_CHECK_EQUAL(tc.size(), 10);

  unsigned int i = 0;
  for (auto track : tc) {
    BOOST_CHECK_EQUAL(i, track.tipIndex());
    track.parameters().setRandom();
    i++;
  }

  BOOST_CHECK_EQUAL(std::distance(tc.begin(), tc.end()), tc.size());
}

BOOST_AUTO_TEST_CASE(IteratorConcept) {
  VectorTrackContainer vtc;
  VectorMultiTrajectory mtj;
  TrackContainer tc{vtc, mtj};

  for (unsigned int i = 0; i < 10; i++) {
    auto t = tc.makeTrack();
    t.tipIndex() = i;
  }
  BOOST_CHECK_EQUAL(tc.size(), 10);
  BOOST_CHECK_EQUAL(std::distance(tc.begin(), tc.end()), tc.size());

  {
    auto it = tc.begin();
    BOOST_CHECK(*it == tc.getTrack(0));
    ++it;
    BOOST_CHECK(*it == tc.getTrack(1));
    it += 1;
    BOOST_CHECK(*it == tc.getTrack(2));
    it -= 1;
    BOOST_CHECK(*it == tc.getTrack(1));
    ++it;
    ++it;
    --it;
    BOOST_CHECK(*it == tc.getTrack(2));
  }
  {
    auto it = tc.begin();
    BOOST_CHECK(*it == tc.getTrack(0));
    std::advance(it, 4);
    BOOST_CHECK(*it == tc.getTrack(4));
    BOOST_CHECK(*(it[-1]) == tc.getTrack(3));
    BOOST_CHECK(*(it[0]) == tc.getTrack(4));
    BOOST_CHECK(*(it[1]) == tc.getTrack(5));
    BOOST_CHECK(*(it - 2) == tc.getTrack(2));
  }

  {
    auto it = tc.begin();
    auto it4 = it + 4;
    auto it5 = it + 5;
    auto it6 = it + 6;

    BOOST_CHECK(it4 < it5);
    BOOST_CHECK(it5 < it6);
    BOOST_CHECK(it4 < it6);

    BOOST_CHECK(it6 > it5);
    BOOST_CHECK(it5 > it4);
    BOOST_CHECK(it6 > it4);

    BOOST_CHECK(it4 <= it4);
    BOOST_CHECK(it4 <= it5);
    BOOST_CHECK(it5 <= it5);
    BOOST_CHECK(it5 <= it6);

    BOOST_CHECK(it6 >= it6);
    BOOST_CHECK(it6 >= it5);
  }
}

BOOST_AUTO_TEST_CASE(ConstCorrectness) {
  VectorTrackContainer vtc;
  VectorMultiTrajectory mtj;
  {
    TrackContainer tc{vtc, mtj};

    for (unsigned int i = 0; i < 10; i++) {
      auto t = tc.makeTrack();
      t.tipIndex() = i;
    }

    unsigned int i = 0;
    for (auto track : tc) {
      BOOST_CHECK_EQUAL(i, track.tipIndex());
      track.parameters().setRandom();
      i++;
    }

    for (const auto track : tc) {
      (void)track;
      // does not compile
      // track.parameters().setRandom();
    }
  }

  ConstVectorTrackContainer cvtc{std::move(vtc)};
  ConstVectorMultiTrajectory cmtj{std::move(mtj)};
  {
    TrackContainer tc{cvtc, cmtj};

    unsigned int i = 0;
    for (auto track : tc) {
      BOOST_CHECK_EQUAL(i, track.tipIndex());
      i++;
      // does not compile
      // track.parameters().setRandom();
    }
  }
}

BOOST_AUTO_TEST_CASE(BuildFromConstRef) {
  VectorTrackContainer mutVtc;
  VectorMultiTrajectory mutMtj;

  TrackContainer mutTc{mutVtc, mutMtj};
  static_assert(!mutTc.ReadOnly, "Unexpectedly read only");

  auto t = mutTc.makeTrack();
  t.appendTrackState();
  t.appendTrackState();
  t.appendTrackState();
  t = mutTc.makeTrack();
  t.appendTrackState();

  BOOST_CHECK_EQUAL(mutTc.size(), 2);
  BOOST_CHECK_EQUAL(mutMtj.size(), 4);

  ConstVectorTrackContainer vtc{std::move(mutVtc)};
  ConstVectorMultiTrajectory mtj{std::move(mutMtj)};

  // moved from
  BOOST_CHECK_EQUAL(mutTc.size(), 0);
  BOOST_CHECK_EQUAL(mutMtj.size(), 0);

  TrackContainer ctc{vtc, mtj};
  static_assert(ctc.ReadOnly, "Unexpectedly not read only");

  // Does not compile:
  // ctc.addTrack();

  BOOST_CHECK_EQUAL(ctc.size(), 2);
  BOOST_CHECK_EQUAL(mtj.size(), 4);

  const auto& cvtc = vtc;
  const auto& cmtj = mtj;

  TrackContainer crtc{cvtc, cmtj};

  BOOST_CHECK_EQUAL(crtc.size(), 2);
  BOOST_CHECK_EQUAL(cmtj.size(), 4);

  // Does not compile: holder deduced to ConstRefHolder, but is not RO
  // const auto& mrvtc = mutVtc;
  // const auto& mrmtj = mutMtj;
  // TrackContainer mrtc{mrvtc, mrmtj};
  // static_assert(ctc.ReadOnly, "Unexpectedly not read only");
  // mrtc.addTrack();
}

BOOST_AUTO_TEST_CASE_TEMPLATE(BuildReadOnly, factory_t, const_holder_types) {
  factory_t factory;
  auto& tc = factory.trackContainer();

  static_assert(std::is_same_v<std::decay_t<decltype(tc)>,
                               typename factory_t::track_container_type>,
                "Incorrect deduction");

  static_assert(std::decay_t<decltype(tc)>::ReadOnly, "Should be read only");
  BOOST_CHECK(tc.ReadOnly);
}

BOOST_AUTO_TEST_CASE_TEMPLATE(DynamicColumns, factory_t, holder_types) {
  factory_t factory;
  auto& tc = factory.trackContainer();

  BOOST_CHECK(!tc.hasColumn("col_a"_hash));
  tc.template addColumn<float>("col_a");
  BOOST_CHECK(tc.hasColumn("col_a"_hash));

  auto t = tc.makeTrack();
  t.template component<float>("col_a") = 5.6f;
  BOOST_CHECK_EQUAL((t.template component<float, "col_a"_hash>()), 5.6f);
}

BOOST_AUTO_TEST_CASE(EnsureDynamicColumns) {
  TrackContainer tc{VectorTrackContainer{}, VectorMultiTrajectory{}};
  tc.addColumn<std::size_t>("counter");
  tc.addColumn<bool>("odd");

  BOOST_CHECK(tc.hasColumn("counter"));
  BOOST_CHECK(tc.hasColumn("odd"));

  TrackContainer tc2{VectorTrackContainer{}, VectorMultiTrajectory{}};

  BOOST_CHECK(!tc2.hasColumn("counter"));
  BOOST_CHECK(!tc2.hasColumn("odd"));

  tc2.ensureDynamicColumns(tc);

  BOOST_CHECK(tc2.hasColumn("counter"));
  BOOST_CHECK(tc2.hasColumn("odd"));
}

BOOST_AUTO_TEST_CASE(AppendTrackState) {
  TrackContainer tc{VectorTrackContainer{}, VectorMultiTrajectory{}};
  auto t = tc.makeTrack();

  std::vector<VectorMultiTrajectory::TrackStateProxy> trackStates;
  trackStates.push_back(t.appendTrackState());
  trackStates.push_back(t.appendTrackState());
  trackStates.push_back(t.appendTrackState());
  trackStates.push_back(t.appendTrackState());
  trackStates.push_back(t.appendTrackState());
  trackStates.push_back(t.appendTrackState());

  BOOST_CHECK_EQUAL(trackStates.size(), t.nTrackStates());

  for (std::size_t i = trackStates.size() - 1; i > 0; i--) {
    BOOST_CHECK_EQUAL(trackStates.at(i).index(), i);
  }
}

BOOST_AUTO_TEST_CASE(ForwardIteration) {
  TrackContainer tc{VectorTrackContainer{}, VectorMultiTrajectory{}};
  {
    // let's create an unrelated track first
    auto t = tc.makeTrack();
    for (std::size_t i = 0; i < 10; i++) {
      t.appendTrackState();
    }
  }

  auto t = tc.makeTrack();

  auto stem = t.appendTrackState();
  t.appendTrackState();
  t.appendTrackState();
  t.appendTrackState();
  t.appendTrackState();

  BOOST_CHECK_THROW(t.trackStates(), std::invalid_argument);
  BOOST_CHECK(!t.innermostTrackState().has_value());

  t.linkForward();

  BOOST_CHECK_EQUAL(t.stemIndex(), stem.index());
  BOOST_CHECK_EQUAL(t.innermostTrackState().value().index(), stem.index());
  t.innermostTrackState()->predicted().setRandom();

  std::vector<IndexType> indices;
  for (const auto& ts : t.trackStatesReversed()) {
    indices.push_back(ts.index());
  }

  std::reverse(indices.begin(), indices.end());

  std::vector<IndexType> act;
  for (auto ts : t.trackStates()) {
    act.push_back(ts.index());
    ts.predicted().setRandom();
  }

  BOOST_CHECK_EQUAL_COLLECTIONS(indices.begin(), indices.end(), act.begin(),
                                act.end());

  t.reverseTrackStates();
  BOOST_CHECK_EQUAL(t.innermostTrackState().value().index(), indices.back());
  t.innermostTrackState()->predicted().setRandom();

  act.clear();
  for (const auto& ts : t.trackStates()) {
    act.push_back(ts.index());
  }

  BOOST_CHECK_EQUAL_COLLECTIONS(indices.rbegin(), indices.rend(), act.begin(),
                                act.end());
}

BOOST_AUTO_TEST_CASE(CalculateQuantities) {
  TrackContainer tc{VectorTrackContainer{}, VectorMultiTrajectory{}};
  auto t = tc.makeTrack();

  auto ts = t.appendTrackState();
  ts.typeFlags().set(MeasurementFlag);

  ts = t.appendTrackState();
  ts.typeFlags().set(OutlierFlag);

  ts = t.appendTrackState();
  ts.typeFlags().set(MeasurementFlag);
  ts.typeFlags().set(SharedHitFlag);

  ts = t.appendTrackState();
  ts.typeFlags().set(HoleFlag);

  ts = t.appendTrackState();
  ts.typeFlags().set(OutlierFlag);

  ts = t.appendTrackState();
  ts.typeFlags().set(HoleFlag);

  ts = t.appendTrackState();
  ts.typeFlags().set(MeasurementFlag);
  ts.typeFlags().set(SharedHitFlag);

  ts = t.appendTrackState();
  ts.typeFlags().set(OutlierFlag);

  calculateTrackQuantities(t);

  BOOST_CHECK_EQUAL(t.nHoles(), 2);
  BOOST_CHECK_EQUAL(t.nMeasurements(), 3);
  BOOST_CHECK_EQUAL(t.nOutliers(), 3);
  BOOST_CHECK_EQUAL(t.nSharedHits(), 2);
}

BOOST_AUTO_TEST_SUITE_END()

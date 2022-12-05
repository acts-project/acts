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

// template <template <typename> class holder_t>
// using track_container_t =
// TrackContainer<VectorTrackContainer, VectorMultiTrajectory, holder_t>;

template <typename track_container_t, typename traj_t,
          template <typename> class holder_t>
struct Factory {};

template <typename track_container_t, typename traj_t>
struct Factory<track_container_t, traj_t, detail_tc::RefHolder> {
  using track_container_type =
      TrackContainer<track_container_t, traj_t, detail_tc::RefHolder>;

  track_container_t vtc;
  traj_t mtj;
  track_container_type tc{vtc, mtj};

  auto& trackContainer() { return tc; }
  auto& trackStateContainer() { return mtj; }
  auto& backend() { return vtc; }
};

template <typename track_container_t, typename traj_t>
struct Factory<track_container_t, traj_t, detail_tc::ValueHolder> {
  using track_container_type =
      TrackContainer<track_container_t, traj_t, detail_tc::ValueHolder>;

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
                                    detail_tc::ValueHolder,
                                    detail_tc::RefHolder, std::shared_ptr>;

using const_holder_types =
    holder_types_t<ConstVectorTrackContainer, ConstVectorMultiTrajectory,
                   detail_tc::ValueHolder, detail_tc::RefHolder,
                   std::shared_ptr>;

}  // namespace

BOOST_AUTO_TEST_SUITE(EventDataTrack)

BOOST_AUTO_TEST_CASE(BuildDefaultHolder) {
  VectorMultiTrajectory mtj{};
  VectorTrackContainer vtc{};
  TrackContainer tc{vtc, mtj};

  static_assert(
      std::is_same_v<decltype(tc),
                     TrackContainer<VectorTrackContainer, VectorMultiTrajectory,
                                    detail_tc::RefHolder>>,
      "Incorrect deduced type");
  BOOST_CHECK_EQUAL(&mtj, &tc.trackStateContainer());
  BOOST_CHECK_EQUAL(&vtc, &tc.container());
  tc.addTrack();
}

BOOST_AUTO_TEST_CASE(BuildValueHolder) {
  {
    VectorMultiTrajectory mtj{};
    VectorTrackContainer vtc{};
    TrackContainer tc{std::move(vtc), std::move(mtj)};
    static_assert(
        std::is_same_v<decltype(tc), TrackContainer<VectorTrackContainer,
                                                    VectorMultiTrajectory,
                                                    detail_tc::ValueHolder>>,
        "Incorrect deduced type");
  }
  {
    TrackContainer tc{VectorTrackContainer{}, VectorMultiTrajectory{}};

    static_assert(
        std::is_same_v<decltype(tc), TrackContainer<VectorTrackContainer,
                                                    VectorMultiTrajectory,
                                                    detail_tc::ValueHolder>>,
        "Incorrect deduced type");
    tc.addTrack();
  }
  {
    VectorMultiTrajectory mtj{};
    VectorTrackContainer vtc{};
    TrackContainer<VectorTrackContainer, VectorMultiTrajectory,
                   detail_tc::ValueHolder>
        tc{vtc, mtj};
    BOOST_CHECK_NE(&mtj, &tc.trackStateContainer());
    BOOST_CHECK_NE(&vtc, &tc.container());
    tc.addTrack();
  }
}

BOOST_AUTO_TEST_CASE(BuildRefHolder) {
  VectorMultiTrajectory mtj{};
  VectorTrackContainer vtc{};
  TrackContainer<VectorTrackContainer, VectorMultiTrajectory,
                 detail_tc::RefHolder>
      tc{vtc, mtj};

  static_assert(
      std::is_same_v<decltype(tc),
                     TrackContainer<VectorTrackContainer, VectorMultiTrajectory,
                                    detail_tc::RefHolder>>,
      "Incorrect deduced type");
  BOOST_CHECK_EQUAL(&mtj, &tc.trackStateContainer());
  BOOST_CHECK_EQUAL(&vtc, &tc.container());
  tc.addTrack();
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

  auto t = tc.getTrack(tc.addTrack());
  t.tipIndex() = ts5.index();

  std::vector<IndexType> act;
  for (auto ts : t.trackStates()) {
    act.push_back(ts.index());
  }

  std::vector<IndexType> exp;
  exp.resize(5);
  std::iota(exp.rbegin(), exp.rend(), 0);
  BOOST_CHECK_EQUAL_COLLECTIONS(act.begin(), act.end(), exp.begin(), exp.end());
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

  auto t = tc.getTrack(tc.addTrack());
  t.template component<float>("col_a") = 5.6f;
  BOOST_CHECK_EQUAL((t.template component<float, "col_a"_hash>()), 5.6f);
}

BOOST_AUTO_TEST_SUITE_END()

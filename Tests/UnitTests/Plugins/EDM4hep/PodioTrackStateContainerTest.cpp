// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/data/test_case.hpp>
#include <boost/test/unit_test.hpp>
#include <boost/test/unit_test_suite.hpp>

#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/EventData/SourceLink.hpp"
#include "Acts/EventData/TrackStatePropMask.hpp"
#include "Acts/EventData/detail/MultiTrajectoryTestsCommon.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/Surfaces/RectangleBounds.hpp"
#include "ActsPlugins/EDM4hep/PodioTrackStateContainer.hpp"
#include "ActsPlugins/EDM4hep/PodioUtil.hpp"
#include "ActsPodioEdm/BoundParametersCollection.h"
#include "ActsPodioEdm/JacobianCollection.h"
#include "ActsPodioEdm/TrackStateCollection.h"

#include <filesystem>

#include <podio/UserDataCollection.h>

namespace {

using namespace Acts;
using namespace ActsPlugins;
using namespace UnitLiterals;
using namespace Acts::detail::Test;
namespace bd = boost::unit_test::data;

std::default_random_engine rng(31415);

class NullHelper : public PodioUtil::ConversionHelper {
 public:
  std::optional<PodioUtil::Identifier> surfaceToIdentifier(
      const Surface& /*surface*/) const override {
    return {};
  }
  const Surface* identifierToSurface(
      PodioUtil::Identifier /*identifier*/) const override {
    return nullptr;
  }

  std::optional<SourceLink> identifierToSourceLink(
      PodioUtil::Identifier /*identifier*/) const override {
    return SourceLink{0};
  }

  std::optional<PodioUtil::Identifier> sourceLinkToIdentifier(
      const SourceLink& /*sourceLink*/) const override {
    return 0;
  }
};

struct MapHelper : public NullHelper {
  std::optional<PodioUtil::Identifier> surfaceToIdentifier(
      const Surface& surface) const override {
    for (auto&& [id, srf] : surfaces) {
      if (srf == &surface) {
        return id;
      }
    }
    return {};
  }
  const Surface* identifierToSurface(PodioUtil::Identifier id) const override {
    auto it = surfaces.find(id);
    if (it == surfaces.end()) {
      return nullptr;
    }

    return it->second;
  }

  std::optional<PodioUtil::Identifier> sourceLinkToIdentifier(
      const SourceLink& sl) const override {
    sourceLinks.push_back(sl);
    return sourceLinks.size() - 1;
  }

  std::optional<SourceLink> identifierToSourceLink(
      PodioUtil::Identifier id) const override {
    return sourceLinks.at(id);
  }

  std::unordered_map<PodioUtil::Identifier, const Surface*> surfaces;
  mutable std::vector<SourceLink> sourceLinks;
};

struct Factory {
  using trajectory_t = MutablePodioTrackStateContainer<>;
  using const_trajectory_t = ConstPodioTrackStateContainer<>;

  MapHelper m_helper;

  MutablePodioTrackStateContainer<> create() {
    return MutablePodioTrackStateContainer{
        m_helper, std::make_unique<ActsPodioEdm::TrackStateCollection>(),
        std::make_unique<ActsPodioEdm::BoundParametersCollection>(),
        std::make_unique<ActsPodioEdm::JacobianCollection>(), nullptr};
  }
};

using CommonTests = MultiTrajectoryTestsCommon<Factory>;

}  // namespace

BOOST_AUTO_TEST_SUITE(PodioTrackStateContainerTest)

BOOST_AUTO_TEST_CASE(Build) {
  CommonTests ct;
  ct.testBuild();
}

BOOST_AUTO_TEST_CASE(ConstCorrectness) {
  // @TODO: Const version can only be non-owning!
}

BOOST_AUTO_TEST_CASE(Clear) {
  CommonTests ct;
  ct.testClear();
}

BOOST_AUTO_TEST_CASE(ApplyWithAbort) {
  CommonTests ct;
  ct.testApplyWithAbort();
}

BOOST_AUTO_TEST_CASE(AddTrackStateWithBitMask) {
  CommonTests ct;
  ct.testAddTrackStateWithBitMask();
}

BOOST_AUTO_TEST_CASE(AddTrackStateComponents) {
  CommonTests ct;
  ct.testAddTrackStateComponents();
}

// assert expected "cross-talk" between trackstate proxies
BOOST_AUTO_TEST_CASE(TrackStateProxyCrossTalk) {
  CommonTests ct;
  ct.testTrackStateProxyCrossTalk(rng);
}

BOOST_AUTO_TEST_CASE(TrackStateReassignment) {
  CommonTests ct;
  ct.testTrackStateReassignment(rng);
}

BOOST_DATA_TEST_CASE(TrackStateProxyStorage, bd::make({1u, 2u}),
                     nMeasurements) {
  CommonTests ct;
  ct.testTrackStateProxyStorage(rng, nMeasurements);
}

BOOST_AUTO_TEST_CASE(TrackStateProxyAllocations) {
  CommonTests ct;
  ct.testTrackStateProxyAllocations(rng);
}

BOOST_AUTO_TEST_CASE(TrackStateProxyGetMask) {
  CommonTests ct;
  ct.testTrackStateProxyGetMask();
}

BOOST_AUTO_TEST_CASE(TrackStateProxyCopy) {
  CommonTests ct;
  ct.testTrackStateProxyCopy(rng);
}

BOOST_AUTO_TEST_CASE(TrackStateCopyDynamicColumns) {
  CommonTests ct;
  ct.testTrackStateCopyDynamicColumns();
}

BOOST_AUTO_TEST_CASE(TrackStateProxyCopyDiffMTJ) {
  CommonTests ct;
  ct.testTrackStateProxyCopyDiffMTJ();
}

BOOST_AUTO_TEST_CASE(ProxyAssignment) {
  CommonTests ct;
  ct.testProxyAssignment();
}

BOOST_AUTO_TEST_CASE(CopyFromConst) {
  CommonTests ct;
  ct.testCopyFromConst();
}

BOOST_AUTO_TEST_CASE(TrackStateProxyShare) {
  CommonTests ct;
  ct.testTrackStateProxyShare(rng);
}

BOOST_AUTO_TEST_CASE(MultiTrajectoryExtraColumns) {
  CommonTests ct;
  ct.testMultiTrajectoryExtraColumns();
}

BOOST_AUTO_TEST_CASE(MultiTrajectoryExtraColumnsRuntime) {
  CommonTests ct;
  ct.testMultiTrajectoryExtraColumnsRuntime();
}

BOOST_AUTO_TEST_CASE(MultiTrajectoryAllocateCalibratedInit) {
  CommonTests ct;
  ct.testMultiTrajectoryAllocateCalibratedInit(rng);
}

BOOST_AUTO_TEST_CASE(WriteToPodioFrame) {
  using namespace HashedStringLiteral;

  MapHelper helper;

  auto tmp_path = std::filesystem::temp_directory_path();
  auto outfile = tmp_path / "trackstates.root";

  BoundVector tv1;
  tv1 << 1, 1, 1, 1, 1, 1;

  BoundVector tv2 = tv1 * 2;
  BoundVector tv3 = tv1 * 3;
  BoundVector tv4 = tv1 * 4;

  BoundMatrix cov1;
  cov1.setOnes();

  BoundMatrix cov2 = cov1 * 2;
  BoundMatrix cov3 = cov1 * 3;
  BoundMatrix cov4 = cov1 * 4;

  auto rBounds = std::make_shared<RectangleBounds>(15, 20);
  auto trf = Transform3::Identity();
  trf.translation().setRandom();
  auto free = Surface::makeShared<PlaneSurface>(trf, rBounds);
  auto reg = Surface::makeShared<PlaneSurface>(trf, rBounds);

  helper.surfaces[666] = reg.get();

  podio::Frame frame;

  MutablePodioTrackStateContainer c{
      helper, std::make_unique<ActsPodioEdm::TrackStateCollection>(),
      std::make_unique<ActsPodioEdm::BoundParametersCollection>(),
      std::make_unique<ActsPodioEdm::JacobianCollection>(), nullptr};
  BOOST_CHECK(!c.hasColumn("int_column"_hash));
  BOOST_CHECK(!c.hasColumn("float_column"_hash));
  c.addColumn<std::int32_t>("int_column");
  c.addColumn<float>("float_column");
  BOOST_CHECK(c.hasColumn("int_column"_hash));
  BOOST_CHECK(c.hasColumn("float_column"_hash));

  {
    auto t1 = c.makeTrackState(TrackStatePropMask::Predicted);
    t1.predicted() = tv1;
    t1.predictedCovariance() = cov1;

    t1.setReferenceSurface(free);

    auto t2 = c.makeTrackState(TrackStatePropMask::All, t1.index());
    t2.predicted() = tv2;
    t2.predictedCovariance() = cov2;

    t2.filtered() = tv3;
    t2.filteredCovariance() = cov3;

    t2.smoothed() = tv4;
    t2.smoothedCovariance() = cov4;

    t2.jacobian() = cov2;

    auto t3 = c.makeTrackState();
    t3.setReferenceSurface(reg);

    t1.component<std::int32_t, "int_column"_hash>() = -11;
    t2.component<std::int32_t, "int_column"_hash>() = 42;
    t3.component<std::int32_t, "int_column"_hash>() = -98;

    t1.component<float, "float_column"_hash>() = -11.2f;
    t2.component<float, "float_column"_hash>() = 42.4f;
    t3.component<float, "float_column"_hash>() = -98.9f;
  }

  c.releaseInto(frame, "test");

  BOOST_CHECK_EQUAL(frame.get("test_trackStates")->size(), 3);
  BOOST_CHECK_EQUAL(frame.get("test_trackStateParameters")->size(), 7);
  BOOST_CHECK_EQUAL(frame.get("test_trackStateJacobians")->size(), 2);
  BOOST_CHECK_NE(frame.get("test_trackStates_extra__int_column"), nullptr);
  BOOST_CHECK_NE(frame.get("test_trackStates_extra__float_column"), nullptr);

  ConstPodioTrackStateContainer cc{helper, frame, "test"};

  BOOST_CHECK_EQUAL(cc.size(), 3);
  BOOST_CHECK(cc.hasColumn("int_column"_hash));
  BOOST_CHECK(cc.hasColumn("float_column"_hash));

  auto t1 = cc.getTrackState(0);
  auto t2 = cc.getTrackState(1);
  auto t3 = cc.getTrackState(2);

  BOOST_CHECK_EQUAL(t2.previous(), 0);

  BOOST_CHECK(t1.hasReferenceSurface());
  BOOST_CHECK(!t2.hasReferenceSurface());
  BOOST_CHECK(t3.hasReferenceSurface());

  auto gctx = GeometryContext::dangerouslyDefaultConstruct();

  const auto& ext = t1.referenceSurface();
  BOOST_CHECK_NE(&ext, free.get());
  BOOST_CHECK_EQUAL(trf.matrix(), ext.localToGlobalTransform(gctx).matrix());
  BOOST_CHECK_EQUAL(free->bounds().type(), ext.bounds().type());
  BOOST_CHECK_EQUAL(free->type(), ext.type());
  const auto* rBounds2 = dynamic_cast<const RectangleBounds*>(&ext.bounds());
  BOOST_REQUIRE_NE(rBounds2, nullptr);
  BOOST_CHECK_EQUAL(rBounds->halfLengthX(), rBounds2->halfLengthX());
  BOOST_CHECK_EQUAL(rBounds->halfLengthY(), rBounds2->halfLengthY());

  BOOST_CHECK_EQUAL(t1.predicted(), tv1);
  BOOST_CHECK_EQUAL(t1.predictedCovariance(), cov1);

  BOOST_CHECK_EQUAL(t2.predicted(), tv2);
  BOOST_CHECK_EQUAL(t2.predictedCovariance(), cov2);
  BOOST_CHECK_EQUAL(t2.filtered(), tv3);
  BOOST_CHECK_EQUAL(t2.filteredCovariance(), cov3);
  BOOST_CHECK_EQUAL(t2.smoothed(), tv4);
  BOOST_CHECK_EQUAL(t2.smoothedCovariance(), cov4);

  BOOST_CHECK_EQUAL(t2.jacobian(), cov2);

  BOOST_CHECK_EQUAL(&t3.referenceSurface(), reg.get());

  BOOST_CHECK_EQUAL((t1.component<std::int32_t, "int_column"_hash>()), -11);
  BOOST_CHECK_EQUAL((t2.component<std::int32_t, "int_column"_hash>()), 42);
  BOOST_CHECK_EQUAL((t3.component<std::int32_t, "int_column"_hash>()), -98);

  BOOST_CHECK_EQUAL((t1.component<float, "float_column"_hash>()), -11.2f);
  BOOST_CHECK_EQUAL((t2.component<float, "float_column"_hash>()), 42.4f);
  BOOST_CHECK_EQUAL((t3.component<float, "float_column"_hash>()), -98.9f);
}

BOOST_AUTO_TEST_CASE(ExternalCollectionSupport) {
  NullHelper helper;

  // Test external collection constructor (non-owning via RefHolder)
  ActsPodioEdm::TrackStateCollection externalTrackStates;
  ActsPodioEdm::BoundParametersCollection externalParams;
  ActsPodioEdm::JacobianCollection externalJacs;

  MutablePodioTrackStateContainer<Acts::RefHolder> externalContainer{
      helper, externalTrackStates, externalParams, externalJacs, nullptr};

  // Add a track state to the external container
  auto ts = externalContainer.addTrackState();
  BOOST_CHECK_EQUAL(ts, 0);
  BOOST_CHECK_EQUAL(externalContainer.size(), 1);
  BOOST_CHECK_EQUAL(externalTrackStates.size(), 1);

  // Test owned collection constructor (default using std::unique_ptr)
  MutablePodioTrackStateContainer ownedContainer{
      helper, std::make_unique<ActsPodioEdm::TrackStateCollection>(),
      std::make_unique<ActsPodioEdm::BoundParametersCollection>(),
      std::make_unique<ActsPodioEdm::JacobianCollection>(), nullptr};
  ts = ownedContainer.addTrackState();
  BOOST_CHECK_EQUAL(ts, 0);
  BOOST_CHECK_EQUAL(ownedContainer.size(), 1);

  // Test that releaseInto works for owned collections
  podio::Frame frame;
  BOOST_CHECK_NO_THROW(ownedContainer.releaseInto(frame, ""));
}

BOOST_AUTO_TEST_CASE(CopyAndMoveConstructors) {
  using namespace HashedStringLiteral;

  MapHelper helper;

  // Create a mutable container with some data
  MutablePodioTrackStateContainer mutableContainer{
      helper, std::make_unique<ActsPodioEdm::TrackStateCollection>(),
      std::make_unique<ActsPodioEdm::BoundParametersCollection>(),
      std::make_unique<ActsPodioEdm::JacobianCollection>(), nullptr};

  // Add some track states
  mutableContainer.addColumn<std::int32_t>("test_column");
  auto ts1 = mutableContainer.addTrackState();
  auto ts2 = mutableContainer.addTrackState();

  BoundVector tv1;
  tv1 << 1, 2, 3, 4, 5, 6;
  mutableContainer.getTrackState(ts1).predicted() = tv1;
  mutableContainer.getTrackState(ts1)
      .template component<std::int32_t, "test_column"_hash>() = 42;

  BoundVector tv2;
  tv2 << 7, 8, 9, 10, 11, 12;
  mutableContainer.getTrackState(ts2).predicted() = tv2;
  mutableContainer.getTrackState(ts2)
      .template component<std::int32_t, "test_column"_hash>() = 99;

  BOOST_CHECK_EQUAL(mutableContainer.size(), 2);

  // Test copy constructor from mutable to const
  ConstPodioTrackStateContainer constContainerCopy{mutableContainer};

  BOOST_CHECK_EQUAL(constContainerCopy.size(), 2);
  BOOST_CHECK(constContainerCopy.hasColumn("test_column"_hash));
  BOOST_CHECK_EQUAL(constContainerCopy.getTrackState(ts1).predicted(), tv1);
  BOOST_CHECK_EQUAL(constContainerCopy.getTrackState(ts2).predicted(), tv2);
  BOOST_CHECK_EQUAL(
      (constContainerCopy.getTrackState(ts1)
           .template component<std::int32_t, "test_column"_hash>()),
      42);
  BOOST_CHECK_EQUAL(
      (constContainerCopy.getTrackState(ts2)
           .template component<std::int32_t, "test_column"_hash>()),
      99);

  // Original mutable container should still be valid
  BOOST_CHECK_EQUAL(mutableContainer.size(), 2);

  // Test move constructor from mutable to const
  // Create a new mutable container for moving
  MutablePodioTrackStateContainer mutableContainer2{
      helper, std::make_unique<ActsPodioEdm::TrackStateCollection>(),
      std::make_unique<ActsPodioEdm::BoundParametersCollection>(),
      std::make_unique<ActsPodioEdm::JacobianCollection>(), nullptr};

  mutableContainer2.addColumn<std::int32_t>("test_column");
  ts1 = mutableContainer2.addTrackState();
  ts2 = mutableContainer2.addTrackState();
  mutableContainer2.getTrackState(ts1).predicted() = tv1;
  mutableContainer2.getTrackState(ts1)
      .template component<std::int32_t, "test_column"_hash>() = 42;
  mutableContainer2.getTrackState(ts2).predicted() = tv2;
  mutableContainer2.getTrackState(ts2)
      .template component<std::int32_t, "test_column"_hash>() = 99;

  // Move construct
  ConstPodioTrackStateContainer constContainerMove{mutableContainer2};

  BOOST_CHECK_EQUAL(constContainerMove.size(), 2);
  BOOST_CHECK(constContainerMove.hasColumn("test_column"_hash));
  BOOST_CHECK_EQUAL(constContainerMove.getTrackState(ts1).predicted(), tv1);
  BOOST_CHECK_EQUAL(constContainerMove.getTrackState(ts2).predicted(), tv2);
  BOOST_CHECK_EQUAL(
      (constContainerMove.getTrackState(ts1)
           .template component<std::int32_t, "test_column"_hash>()),
      42);
  BOOST_CHECK_EQUAL(
      (constContainerMove.getTrackState(ts2)
           .template component<std::int32_t, "test_column"_hash>()),
      99);
}

// Tests for the link-based uncalibrated source link storage mode, which is
// activated by passing a non-null TrackStateHitLinkCollection* to the
// container constructor.

BOOST_AUTO_TEST_CASE(UncalibratedSourceLinkLinkModeBasic) {
  NullHelper helper;

  ActsPodioEdm::TrackStateCollection trackStates;
  ActsPodioEdm::BoundParametersCollection params;
  ActsPodioEdm::JacobianCollection jacs;
  ActsPodioEdm::TrackerHitLocalCollection hitsCollection;
  ActsPodioEdm::TrackStateHitLinkCollection linksCollection;

  // Populate the hits collection with two distinct hits
  auto hit1 = hitsCollection.create();
  hit1.setCellID(12345);
  auto hit2 = hitsCollection.create();
  hit2.setCellID(67890);

  MutablePodioTrackStateContainer<Acts::RefHolder> tsc{
      helper, trackStates, params, jacs, &linksCollection};

  auto i0 = tsc.addTrackState();
  auto i1 = tsc.addTrackState();
  auto i2 = tsc.addTrackState();  // this one gets no source link

  // Assign source links to the first two track states.
  tsc.getTrackState(i0).setUncalibratedSourceLink(
      Acts::SourceLink{hitsCollection.at(0)});
  tsc.getTrackState(i1).setUncalibratedSourceLink(
      Acts::SourceLink{hitsCollection.at(1)});

  // hasUncalibratedSourceLink must reflect link presence
  BOOST_CHECK(tsc.getTrackState(i0).hasUncalibratedSourceLink());
  BOOST_CHECK(tsc.getTrackState(i1).hasUncalibratedSourceLink());
  BOOST_CHECK(!tsc.getTrackState(i2).hasUncalibratedSourceLink());

  // getUncalibratedSourceLink must return the correct TrackerHitLocal
  {
    auto sl = tsc.getTrackState(i0).getUncalibratedSourceLink();
    const auto* hit = sl.getPtr<ActsPodioEdm::TrackerHitLocal>();
    BOOST_REQUIRE_NE(hit, nullptr);
    BOOST_CHECK_EQUAL(hit->getCellID(), 12345u);
  }
  {
    auto sl = tsc.getTrackState(i1).getUncalibratedSourceLink();
    const auto* hit = sl.getPtr<ActsPodioEdm::TrackerHitLocal>();
    BOOST_REQUIRE_NE(hit, nullptr);
    BOOST_CHECK_EQUAL(hit->getCellID(), 67890u);
  }

  // Requesting a source link for a state with no link must throw
  BOOST_CHECK_THROW(tsc.getTrackState(i2).getUncalibratedSourceLink(),
                    std::invalid_argument);

  // Updating an existing link: point i0 to hit2 instead of hit1
  tsc.getTrackState(i0).setUncalibratedSourceLink(
      Acts::SourceLink{hitsCollection.at(1)});
  {
    auto sl = tsc.getTrackState(i0).getUncalibratedSourceLink();
    const auto* hit = sl.getPtr<ActsPodioEdm::TrackerHitLocal>();
    BOOST_REQUIRE_NE(hit, nullptr);
    BOOST_CHECK_EQUAL(hit->getCellID(), 67890u);
  }

  // The link count must not grow when updating — i0 was updated, not duplicated
  BOOST_CHECK_EQUAL(linksCollection.size(), 2u);
}

BOOST_AUTO_TEST_CASE(InvalidHitsLinksCombination) {
  NullHelper helper;

  ActsPodioEdm::TrackStateHitLinkCollection linksCollection;

  // Constructor validity: only the (optional) TrackStateHitLinkCollection
  // pointer controls whether the container uses link-based or identifier-based
  // uncalibrated source-link storage.
  BOOST_CHECK_NO_THROW((MutablePodioTrackStateContainer{
      helper, std::make_unique<ActsPodioEdm::TrackStateCollection>(),
      std::make_unique<ActsPodioEdm::BoundParametersCollection>(),
      std::make_unique<ActsPodioEdm::JacobianCollection>(), nullptr}));

  BOOST_CHECK_NO_THROW((MutablePodioTrackStateContainer{
      helper, std::make_unique<ActsPodioEdm::TrackStateCollection>(),
      std::make_unique<ActsPodioEdm::BoundParametersCollection>(),
      std::make_unique<ActsPodioEdm::JacobianCollection>(), &linksCollection}));
}

BOOST_AUTO_TEST_CASE(UncalibratedSourceLinkLinkModeRoundTrip) {
  // Verify that ConstPodioTrackStateContainer constructed from the same
  // in-memory collections as its mutable counterpart can correctly navigate the
  // links that were set by the mutable container.  The copy-from-mutable
  // constructor is used so that both containers share the same underlying
  // collections.
  NullHelper helper;

  ActsPodioEdm::TrackStateCollection trackStates;
  ActsPodioEdm::BoundParametersCollection params;
  ActsPodioEdm::JacobianCollection jacs;
  ActsPodioEdm::TrackerHitLocalCollection hitsCollection;
  ActsPodioEdm::TrackStateHitLinkCollection linksCollection;

  auto hit1 = hitsCollection.create();
  hit1.setCellID(11111);
  auto hit2 = hitsCollection.create();
  hit2.setCellID(22222);

  MutablePodioTrackStateContainer tsc{
      helper, trackStates, params, jacs, &linksCollection};

  auto i0 = tsc.addTrackState();
  auto i1 = tsc.addTrackState();
  tsc.addTrackState();  // third state intentionally has no source link

  tsc.getTrackState(i0).setUncalibratedSourceLink(
      Acts::SourceLink{hitsCollection.at(0)});
  tsc.getTrackState(i1).setUncalibratedSourceLink(
      Acts::SourceLink{hitsCollection.at(1)});

  // Construct a const container from the mutable one.  The copy constructor
  // passes through m_hits and m_links so the const container is in link-based
  // mode and shares the same underlying collections.
  ConstPodioTrackStateContainer cc{tsc};

  BOOST_CHECK_EQUAL(cc.size(), 3u);

  BOOST_CHECK(cc.getTrackState(0).hasUncalibratedSourceLink());
  BOOST_CHECK(cc.getTrackState(1).hasUncalibratedSourceLink());
  BOOST_CHECK(!cc.getTrackState(2).hasUncalibratedSourceLink());

  {
    auto sl = cc.getTrackState(0).getUncalibratedSourceLink();
    const auto* hit = sl.getPtr<ActsPodioEdm::TrackerHitLocal>();
    BOOST_REQUIRE_NE(hit, nullptr);
    BOOST_CHECK_EQUAL(hit->getCellID(), 11111u);
  }
  {
    auto sl = cc.getTrackState(1).getUncalibratedSourceLink();
    const auto* hit = sl.getPtr<ActsPodioEdm::TrackerHitLocal>();
    BOOST_REQUIRE_NE(hit, nullptr);
    BOOST_CHECK_EQUAL(hit->getCellID(), 22222u);
  }

  // Requesting a source link for the unlinked state must throw
  BOOST_CHECK_THROW(cc.getTrackState(2).getUncalibratedSourceLink(),
                    std::invalid_argument);
}

BOOST_AUTO_TEST_SUITE_END()

// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/test/data/test_case.hpp>
#include <boost/test/tools/old/interface.hpp>
#include <boost/test/unit_test.hpp>
#include <boost/test/unit_test_suite.hpp>

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/Definitions/Units.hpp"
#include "Acts/EventData/MultiTrajectory.hpp"
#include "Acts/EventData/VectorMultiTrajectory.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Plugins/Podio/PodioTrackContainer.hpp"
#include "Acts/Plugins/Podio/PodioTrackStateContainer.hpp"
#include "Acts/Plugins/Podio/PodioUtil.hpp"
#include "Acts/Surfaces/AnnulusBounds.hpp"
#include "Acts/Surfaces/ConeSurface.hpp"
#include "Acts/Surfaces/ConvexPolygonBounds.hpp"
#include "Acts/Surfaces/CylinderBounds.hpp"
#include "Acts/Surfaces/CylinderSurface.hpp"
#include "Acts/Surfaces/DiamondBounds.hpp"
#include "Acts/Surfaces/DiscBounds.hpp"
#include "Acts/Surfaces/DiscSurface.hpp"
#include "Acts/Surfaces/DiscTrapezoidBounds.hpp"
#include "Acts/Surfaces/EllipseBounds.hpp"
#include "Acts/Surfaces/PerigeeSurface.hpp"
#include "Acts/Surfaces/PlanarBounds.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/Surfaces/RadialBounds.hpp"
#include "Acts/Surfaces/RectangleBounds.hpp"
#include "Acts/Surfaces/StrawSurface.hpp"
#include "Acts/Surfaces/SurfaceBounds.hpp"
#include "Acts/Utilities/Helpers.hpp"
#include "ActsPodioEdm/Surface.h"
#include <ActsPodioEdm/TrackCollection.h>

#include <algorithm>
#include <iterator>
#include <memory>
#include <random>
#include <stdexcept>

using namespace Acts;
using namespace Acts::UnitLiterals;
using namespace Acts::HashedStringLiteral;
BOOST_AUTO_TEST_SUITE(PodioTrackConversion)

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

  SourceLink identifierToSourceLink(
      PodioUtil::Identifier /*identifier*/) const override {
    return SourceLink{0};
  }

  PodioUtil::Identifier sourceLinkToIdentifier(
      const SourceLink& /*sourceLink*/) override {
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

  std::unordered_map<PodioUtil::Identifier, const Surface*> surfaces;
};

BOOST_AUTO_TEST_CASE(ConvertSurface) {
  auto rBounds = std::make_shared<RectangleBounds>(15, 20);

  auto trf = Transform3::Identity();
  trf.translation().setRandom();

  auto free = Acts::Surface::makeShared<PlaneSurface>(trf, rBounds);

  NullHelper helper;
  auto surface = PodioUtil::convertSurfaceToPodio(helper, *free);

  auto free2 = PodioUtil::convertSurfaceFromPodio(helper, surface);

  Acts::GeometryContext gctx;

  BOOST_REQUIRE(free2);
  BOOST_CHECK_EQUAL(free->type(), free2->type());
  BOOST_CHECK_EQUAL(free->bounds().type(), free2->bounds().type());
  BOOST_CHECK_EQUAL(free->center(gctx), free2->center(gctx));

  const auto* rBounds2 = dynamic_cast<const RectangleBounds*>(&free2->bounds());
  BOOST_REQUIRE_NE(rBounds2, nullptr);

  BOOST_CHECK_EQUAL(rBounds2->halfLengthX(), rBounds->halfLengthX());
  BOOST_CHECK_EQUAL(rBounds2->halfLengthY(), rBounds->halfLengthY());

  // this could probably use some more complete checks
}

BOOST_AUTO_TEST_CASE(ConvertTrack) {
  auto rBounds = std::make_shared<RectangleBounds>(15, 20);
  auto trf = Transform3::Identity();
  trf.translation().setRandom();
  auto free = Acts::Surface::makeShared<PlaneSurface>(trf, rBounds);

  MapHelper helper;

  auto refCov = BoundMatrix::Random().eval();

  podio::Frame frame;

  {
    Acts::MutablePodioTrackStateContainer tsc{helper};
    Acts::MutablePodioTrackContainer ptc{helper};
    ActsPodioEdm::TrackCollection& tracks = ptc.trackCollection();

    Acts::TrackContainer tc{ptc, tsc};

    BOOST_CHECK(!tc.hasColumn("int_column"_hash));
    BOOST_CHECK(!tc.hasColumn("float_column"_hash));
    tc.addColumn<int32_t>("int_column");
    tc.addColumn<float>("float_column");
    BOOST_CHECK(tc.hasColumn("int_column"_hash));
    BOOST_CHECK(tc.hasColumn("float_column"_hash));

    BOOST_CHECK_EQUAL(tc.size(), 0);

    auto t = tc.getTrack(tc.addTrack());
    BOOST_CHECK_EQUAL(t.tipIndex(), MultiTrajectoryTraits::kInvalid);

    BOOST_CHECK_EQUAL(tsc.size(), 0);
    auto ts1 = t.appendTrackState();
    auto ts2 = t.appendTrackState();
    auto ts3 = t.appendTrackState();
    BOOST_CHECK_EQUAL(tsc.size(), 3);
    BOOST_CHECK_EQUAL(ts1.index(), 0);
    BOOST_CHECK_EQUAL(ts2.index(), 1);
    BOOST_CHECK_EQUAL(ts3.index(), 2);

    BOOST_CHECK_EQUAL(t.nTrackStates(), 3);
    BOOST_CHECK_EQUAL(t.tipIndex(), 2);

    BOOST_CHECK_EQUAL(tc.size(), 1);

    auto pTrack = tracks.at(0);
    BOOST_CHECK_EQUAL(pTrack.data().tipIndex, 2);

    t.parameters() << 1, 2, 3, 4, 5, 6;
    Eigen::Map<BoundVector> pars{pTrack.data().parameters.data()};
    BoundVector bv;
    bv << 1, 2, 3, 4, 5, 6;
    BOOST_CHECK_EQUAL(pars, bv);

    t.covariance() = refCov;

    Eigen::Map<const BoundMatrix> cov{pTrack.data().covariance.data()};
    BOOST_CHECK_EQUAL(refCov, cov);

    t.nMeasurements() = 17;
    BOOST_CHECK_EQUAL(pTrack.data().nMeasurements, 17);

    t.nHoles() = 34;
    BOOST_CHECK_EQUAL(pTrack.data().nHoles, 34);

    t.chi2() = 882.3f;
    BOOST_CHECK_EQUAL(pTrack.data().chi2, 882.3f);

    t.nDoF() = 9;
    BOOST_CHECK_EQUAL(pTrack.data().ndf, 9);

    t.nOutliers() = 77;
    BOOST_CHECK_EQUAL(pTrack.data().nOutliers, 77);

    t.nSharedHits() = 99;
    BOOST_CHECK_EQUAL(pTrack.data().nSharedHits, 99);

    Acts::GeometryContext gctx;
    t.setReferenceSurface(free);
    const auto& free2 = t.referenceSurface();
    BOOST_CHECK_EQUAL(free->center(gctx), free2.center(gctx));

    const auto* rBounds2 =
        dynamic_cast<const RectangleBounds*>(&free2.bounds());
    BOOST_REQUIRE_NE(rBounds2, nullptr);

    BOOST_CHECK_EQUAL(rBounds2->halfLengthX(), rBounds->halfLengthX());
    BOOST_CHECK_EQUAL(rBounds2->halfLengthY(), rBounds->halfLengthY());

    BOOST_CHECK_EQUAL(pTrack.getReferenceSurface().identifier,
                      PodioUtil::kNoIdentifier);

    auto t2 = tc.getTrack(tc.addTrack());
    auto t3 = tc.getTrack(tc.addTrack());
    BOOST_CHECK_EQUAL(tc.size(), 3);

    // Register surface "with the detector"
    helper.surfaces[666] = free.get();
    t2.setReferenceSurface(free);
    auto pTrack2 = tracks.at(1);
    BOOST_CHECK_EQUAL(pTrack2.getReferenceSurface().identifier, 666);

    t.component<int32_t, "int_column"_hash>() = -11;
    t2.component<int32_t, "int_column"_hash>() = 42;
    t3.component<int32_t, "int_column"_hash>() = -98;

    t.component<float, "float_column"_hash>() = -11.2f;
    t2.component<float, "float_column"_hash>() = 42.4f;
    t3.component<float, "float_column"_hash>() = -98.9f;

    ptc.releaseInto(frame);
    tsc.releaseInto(frame);

    BOOST_REQUIRE_NE(frame.get("tracks"), nullptr);
    BOOST_CHECK_EQUAL(frame.get("tracks")->size(), 3);
    BOOST_REQUIRE_NE(frame.get("tracks_extra__int_column"), nullptr);
    BOOST_REQUIRE_NE(frame.get("tracks_extra__float_column"), nullptr);

    BOOST_REQUIRE_NE(frame.get("trackStates"), nullptr);
    BOOST_CHECK_EQUAL(frame.get("trackStates")->size(), 3);
  }

  {
    Acts::ConstPodioTrackStateContainer tsc{helper, frame};
    Acts::ConstPodioTrackContainer ptc{helper, frame};
    // const ActsPodioEdm::TrackCollection& tracks = ptc.trackCollection();

    Acts::TrackContainer tc{ptc, tsc};

    BOOST_CHECK(tc.hasColumn("int_column"_hash));
    BOOST_CHECK(tc.hasColumn("float_column"_hash));

    BOOST_CHECK_EQUAL(tc.size(), 3);

    auto t = tc.getTrack(0);
    const auto& freeRecreated = t.referenceSurface();
    // Not the exact same surface, it's recreated from values
    BOOST_CHECK_NE(free.get(), &freeRecreated);

    BOOST_CHECK_EQUAL(t.nMeasurements(), 17);

    BOOST_CHECK_EQUAL(t.nHoles(), 34);

    BOOST_CHECK_EQUAL(t.chi2(), 882.3f);

    BOOST_CHECK_EQUAL(t.nDoF(), 9);

    BOOST_CHECK_EQUAL(t.nOutliers(), 77);

    BOOST_CHECK_EQUAL(t.nSharedHits(), 99);

    BOOST_CHECK_EQUAL(t.tipIndex(), 2);
    BOOST_CHECK_EQUAL(t.nTrackStates(), 3);

    auto t2 = tc.getTrack(1);
    // Is the exact same surface, because it's looked up in the "detector"
    BOOST_CHECK_EQUAL(free.get(), &t2.referenceSurface());
    BoundVector bv;
    bv << 1, 2, 3, 4, 5, 6;
    BOOST_CHECK_EQUAL(t.parameters(), bv);

    BOOST_CHECK_EQUAL(t.covariance(), refCov);

    auto t3 = tc.getTrack(2);
    BOOST_CHECK(!t3.hasReferenceSurface());

    BOOST_CHECK_EQUAL((t.component<int32_t, "int_column"_hash>()), -11);
    BOOST_CHECK_EQUAL((t2.component<int32_t, "int_column"_hash>()), 42);
    BOOST_CHECK_EQUAL((t3.component<int32_t, "int_column"_hash>()), -98);

    BOOST_CHECK_EQUAL((t.component<float, "float_column"_hash>()), -11.2f);
    BOOST_CHECK_EQUAL((t2.component<float, "float_column"_hash>()), 42.4f);
    BOOST_CHECK_EQUAL((t3.component<float, "float_column"_hash>()), -98.9f);
  }
}

// @TODO: Add ensure dynamic columns

BOOST_AUTO_TEST_SUITE_END()

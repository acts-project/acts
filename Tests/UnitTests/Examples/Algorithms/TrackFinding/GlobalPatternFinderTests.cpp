// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/Units.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/Surfaces/RectangleBounds.hpp"
#include "Acts/Utilities/Helpers.hpp"
#include "ActsExamples/EventData/MuonSpacePoint.hpp"
#include "ActsExamples/TrackFinding/GlobalPatternFinderAlgorithm.hpp"

#include <cstdint>
#include <memory>
#include <numbers>

using namespace Acts::UnitLiterals;

namespace ActsTests {

namespace {

using MuonId = ActsExamples::MuonSpacePoint::MuonId;
using ActsExamples::ExpandedSector;
using ActsExamples::HitPayload;
using ActsExamples::LayerIndex;
using ActsExamples::PatternTopology;
using ActsExamples::SeedSelector;
using ActsExamples::StIndex;

/// @brief Surface every synthetic hit is placed on. The pattern finder takes the
///        local measurement axes from it.
std::shared_ptr<Acts::Surface> testSurface() {
  return Acts::Surface::makeShared<Acts::PlaneSurface>(
      Acts::Transform3::Identity(),
      std::make_shared<Acts::RectangleBounds>(1._m, 1._m));
}

/// @brief Build a space point with the given chamber & measurement flags
/// @param station: Station name of the chamber
/// @param tech: Technology of the measurement
/// @param layer: Layer number in the sector frame
/// @param measEta: Whether the space point measures the precision coordinate
/// @param measPhi: Whether the space point measures the non-precision coordinate
/// @param pos: Position of the space point in the sector frame
ActsExamples::MuonSpacePoint makeSpacePoint(MuonId::StationName station,
                                            MuonId::TechField tech,
                                            std::uint8_t layer, bool measEta,
                                            bool measPhi,
                                            const Acts::Vector3& pos) {
  MuonId id{};
  id.setChamber(station, MuonId::DetSide::A, 1u, tech);
  id.setLayAndCh(layer, 1u);
  id.setCoordFlags(measEta, measPhi, false);

  ActsExamples::MuonSpacePoint sp{};
  sp.setId(id);
  sp.setGeometryId(Acts::GeometryIdentifier{layer + 1u});
  sp.defineCoordinates(Acts::Vector3{pos}, Acts::Vector3::UnitX(),
                       Acts::Vector3::UnitY());
  // The covariance is ordered as precision, non-precision & time
  sp.setCovariance(0.09, 0., 0.);
  return sp;
}

}  // namespace

BOOST_AUTO_TEST_SUITE(GlobalPatternFinderSuite)

BOOST_AUTO_TEST_CASE(expanded_sector) {
  using SectorProjector = ExpandedSector::SectorProjector;

  // The expanded sector is twice the ms sector plus the projector
  const ExpandedSector centreOfOne{1u, SectorProjector::center};
  BOOST_CHECK_EQUAL(centreOfOne.sector(), 2);
  BOOST_CHECK_EQUAL(centreOfOne.msSector(), 1u);
  BOOST_CHECK(centreOfOne.projector() == SectorProjector::center);
  BOOST_CHECK_EQUAL(centreOfOne.adjacentMsSector(), 1u);
  // Sector 1 is centred at phi = 0
  BOOST_CHECK_SMALL(centreOfOne.phi(), 1.e-9);
  BOOST_CHECK(centreOfOne.insideSector(0.));

  // A hit at phi = 0 is well inside sector 1 and not in an overlap
  const ExpandedSector fromPhi{0.};
  BOOST_CHECK_EQUAL(fromPhi.sector(), centreOfOne.sector());

  // Sector 3 is centred at pi / 4
  const ExpandedSector centreOfThree{3u, SectorProjector::center};
  BOOST_CHECK_CLOSE(centreOfThree.phi(), std::numbers::pi / 4., 1.e-6);
  BOOST_CHECK(!centreOfThree.insideSector(0.));

  // Neighbouring expanded sectors differ by one unit
  BOOST_CHECK(centreOfOne.isNeighbour(
      ExpandedSector{static_cast<ExpandedSector::Index_t>(3)}));
  BOOST_CHECK(!centreOfOne.isNeighbour(centreOfThree));

  // The expanded-sector index wraps at nExpanded (= 32). Sectors 0 and 31
  // meet at the MS boundary and must count as neighbours.
  const ExpandedSector wrapLo{static_cast<ExpandedSector::Index_t>(0)};
  const ExpandedSector wrapHi{static_cast<ExpandedSector::Index_t>(31)};
  BOOST_CHECK(wrapLo.isNeighbour(wrapHi));
  BOOST_CHECK(wrapHi.isNeighbour(wrapLo));
  BOOST_CHECK(!wrapLo.isNeighbour(
      ExpandedSector{static_cast<ExpandedSector::Index_t>(2)}));

  // The overlap between two sectors is expressed as the right overlap of the
  // lower one, which is why the left overlap maps onto the sector below
  const ExpandedSector leftOfTwo{2u, SectorProjector::leftOverlap};
  BOOST_CHECK_EQUAL(leftOfTwo.msSector(), 1u);
  BOOST_CHECK(leftOfTwo.projector() == SectorProjector::rightOverlap);
  BOOST_CHECK_EQUAL(leftOfTwo.adjacentMsSector(), 2u);
  // An overlap is narrower than a sector
  BOOST_CHECK_LT(leftOfTwo.sectorSize(), centreOfOne.sectorSize());
}

BOOST_AUTO_TEST_CASE(station_helpers) {
  using ActsExamples::isBarrel;
  using ActsExamples::isPrecisionHit;
  using ActsExamples::toLayerIndex;
  using ActsExamples::toStationIndex;
  using enum MuonId::StationName;

  BOOST_CHECK(toStationIndex(BIL) == StIndex::BI);
  BOOST_CHECK(toStationIndex(BIS) == StIndex::BI);
  BOOST_CHECK(toStationIndex(BML) == StIndex::BM);
  BOOST_CHECK(toStationIndex(BEE) == StIndex::BE);
  BOOST_CHECK(toStationIndex(EMS) == StIndex::EM);
  BOOST_CHECK(toStationIndex(EEL) == StIndex::EE);

  BOOST_CHECK(toLayerIndex(StIndex::BI) == LayerIndex::Inner);
  BOOST_CHECK(toLayerIndex(StIndex::EM) == LayerIndex::Middle);
  BOOST_CHECK(toLayerIndex(StIndex::BO) == LayerIndex::Outer);
  BOOST_CHECK(toLayerIndex(StIndex::BE) == LayerIndex::BarrelExtended);
  BOOST_CHECK(toLayerIndex(StIndex::EE) == LayerIndex::Extended);

  BOOST_CHECK(isBarrel(StIndex::BM));
  BOOST_CHECK(!isBarrel(StIndex::EM));

  // Mdt & micromega hits are precision hits
  BOOST_CHECK(isPrecisionHit(makeSpacePoint(
      BIL, MuonId::TechField::Mdt, 1u, true, false, Acts::Vector3::Zero())));
  BOOST_CHECK(isPrecisionHit(makeSpacePoint(
      EIL, MuonId::TechField::Mm, 1u, true, false, Acts::Vector3::Zero())));
  // Rpc / Tgc strips are trigger hits
  BOOST_CHECK(!isPrecisionHit(makeSpacePoint(
      BML, MuonId::TechField::Rpc, 1u, true, false, Acts::Vector3::Zero())));
  // An sTgc hit counts as precision only if it measures eta alone. The channel
  // type is not exported, so a combined hit cannot be told from a pad
  BOOST_CHECK(isPrecisionHit(makeSpacePoint(
      EIL, MuonId::TechField::sTgc, 1u, true, false, Acts::Vector3::Zero())));
  BOOST_CHECK(!isPrecisionHit(makeSpacePoint(
      EIL, MuonId::TechField::sTgc, 1u, true, true, Acts::Vector3::Zero())));
}

BOOST_AUTO_TEST_CASE(hit_payload_and_topology) {
  const auto gctx = Acts::GeometryContext::dangerouslyDefaultConstruct();
  const auto surface = testSurface();

  ActsExamples::MuonSpacePointBucket bucket{};
  bucket.push_back(makeSpacePoint(MuonId::StationName::BIL,
                                  MuonId::TechField::Mdt, 1u, true, false,
                                  Acts::Vector3{0., 10., 100.}));
  bucket.push_back(makeSpacePoint(MuonId::StationName::BIL,
                                  MuonId::TechField::Mdt, 2u, true, false,
                                  Acts::Vector3{0., 20., 130.}));
  bucket.push_back(makeSpacePoint(MuonId::StationName::BIL,
                                  MuonId::TechField::Mdt, 2u, true, false,
                                  Acts::Vector3{0., 25., 130.}));

  // A unit transform keeps the sector frame, so the global position is the
  // position the space point was built with
  const Acts::Transform3 toGlobal{Acts::Transform3::Identity()};
  const HitPayload first{gctx, &bucket[0], &bucket, toGlobal, surface.get()};
  const HitPayload second{gctx, &bucket[1], &bucket, toGlobal, surface.get()};
  const HitPayload sameLayerAsSecond{gctx, &bucket[2], &bucket, toGlobal,
                                     surface.get()};

  BOOST_CHECK(first.globalPosition(gctx).isApprox(bucket[0].localPosition()));
  BOOST_CHECK(first.spacePoint() == &bucket[0]);
  BOOST_CHECK(first.station == StIndex::BI);
  BOOST_CHECK_EQUAL(static_cast<unsigned>(first.locLayer),
                    static_cast<unsigned>(bucket[0].id().detLayer()));
  BOOST_CHECK(first.isPrecision());
  BOOST_CHECK(first == first);
  BOOST_CHECK(!(first == second));

  // The hits are grouped by station and ordered by layer
  BOOST_CHECK_EQUAL(static_cast<int>(PatternTopology::groupIndex(first)),
                    static_cast<int>(Acts::toUnderlying(StIndex::BI)));
  BOOST_CHECK(PatternTopology::layerSorter(first, second));
  BOOST_CHECK(!PatternTopology::layerSorter(second, first));
  BOOST_CHECK(!PatternTopology::sameLayer(first, second));
  BOOST_CHECK(PatternTopology::sameLayer(second, sameLayerAsSecond));
  // Two hits on the same layer are ordered by their position along the layer
  BOOST_CHECK(PatternTopology::layerSorter(second, sameLayerAsSecond));
}

BOOST_AUTO_TEST_CASE(seed_selector) {
  const auto gctx = Acts::GeometryContext::dangerouslyDefaultConstruct();
  const auto surface = testSurface();
  const Acts::Transform3 toGlobal{Acts::Transform3::Identity()};

  ActsExamples::MuonSpacePointBucket bucket{};
  bucket.push_back(makeSpacePoint(MuonId::StationName::BML,
                                  MuonId::TechField::Rpc, 1u, true, false,
                                  Acts::Vector3{0., 10., 100.}));
  bucket.push_back(makeSpacePoint(MuonId::StationName::BML,
                                  MuonId::TechField::Mdt, 2u, true, false,
                                  Acts::Vector3{0., 20., 130.}));
  bucket.push_back(makeSpacePoint(MuonId::StationName::BIL,
                                  MuonId::TechField::Rpc, 3u, true, false,
                                  Acts::Vector3{0., 30., 160.}));
  bucket.push_back(makeSpacePoint(MuonId::StationName::BOL,
                                  MuonId::TechField::Rpc, 4u, true, false,
                                  Acts::Vector3{0., 40., 190.}));

  const HitPayload middleStrip{gctx, &bucket[0], &bucket, toGlobal,
                               surface.get()};
  const HitPayload middleStraw{gctx, &bucket[1], &bucket, toGlobal,
                               surface.get()};
  const HitPayload innerStrip{gctx, &bucket[2], &bucket, toGlobal,
                              surface.get()};
  const HitPayload outerStrip{gctx, &bucket[3], &bucket, toGlobal,
                              surface.get()};

  SeedSelector::Config cfg{};
  const SeedSelector selector{SeedSelector::Config{cfg}};

  // Seeds are taken from the middle & outer station, Mdt hits are excluded
  BOOST_CHECK(selector.goodForSeeding(middleStrip));
  BOOST_CHECK(selector.goodForSeeding(outerStrip));
  BOOST_CHECK(!selector.goodForSeeding(middleStraw));
  BOOST_CHECK(!selector.goodForSeeding(innerStrip));

  // A seed's windows are the row of its station in the table of its |eta| bin
  const auto row = [&](const HitPayload& seed) {
    return selector.m_thetaWindows[Acts::toUnderlying(selector.etaBin(seed))]
                                  [Acts::toUnderlying(seed.station)];
  };
  BOOST_CHECK(selector.thetaSearchWindow(middleStrip) == row(middleStrip));
  BOOST_CHECK(selector.thetaSearchWindow(outerStrip) == row(outerStrip));

  // The tables are symmetric, and a barrel middle seed sees the neighbouring
  // barrel stations but not the endcap
  using enum StIndex;
  constexpr auto idx = [](StIndex st) { return Acts::toUnderlying(st); };
  for (const SeedSelector::Matrix_t& m : selector.m_thetaWindows) {
    for (std::size_t a = 0; a < m.size(); ++a) {
      for (std::size_t b = 0; b < m.size(); ++b) {
        BOOST_CHECK_EQUAL(m[a][b], m[b][a]);
      }
    }
  }
  const SeedSelector::Row_t& barrelMiddle{
      selector.m_thetaWindows[Acts::toUnderlying(SeedSelector::EtaBin::Barrel)]
                             [idx(BM)]};
  BOOST_CHECK_EQUAL(barrelMiddle[idx(BM)], cfg.barrelSameStation);
  BOOST_CHECK_EQUAL(barrelMiddle[idx(BI)], cfg.barrelAdjacent);
  BOOST_CHECK_EQUAL(barrelMiddle[idx(BO)], cfg.barrelAdjacent);
  BOOST_CHECK_EQUAL(barrelMiddle[idx(EM)], 0.);
  BOOST_CHECK_CLOSE(selector.maxThetaWindow(), 2.6_degree, 1.e-6);

  // Both toggles open up the corresponding seeds
  SeedSelector::Config mdtCfg{cfg};
  mdtCfg.seedFromMdt = true;
  mdtCfg.seedFromInner = true;
  const SeedSelector fromEverything{SeedSelector::Config{mdtCfg}};
  BOOST_CHECK(fromEverything.goodForSeeding(middleStraw));
  BOOST_CHECK(fromEverything.goodForSeeding(innerStrip));
  BOOST_CHECK(fromEverything.thetaSearchWindow(innerStrip) ==
              fromEverything.m_thetaWindows[Acts::toUnderlying(
                  fromEverything.etaBin(innerStrip))][idx(BI)]);
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace ActsTests

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
#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/EventData/MuonHoughMaximum.hpp"
#include "ActsExamples/EventData/MuonSegment.hpp"
#include "ActsExamples/EventData/MuonSpacePoint.hpp"
#include "ActsExamples/Framework/AlgorithmContext.hpp"
#include "ActsExamples/Framework/DataHandle.hpp"
#include "ActsExamples/Framework/WhiteBoard.hpp"
#include "ActsExamples/TrackFinding/MuonHoughSeeder.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <string>
#include <utility>
#include <vector>

namespace ActsTests {

namespace {

struct DriftCircle {
  double y = 0.0;
  double z = 0.0;
  double rDrift = 0.0;
  double rDriftError = 0.0;
};

/// @brief Prepare MuonId
ActsExamples::MuonSpacePoint::MuonId makeMdtEtaId(std::uint8_t layer,
                                                  std::uint16_t channel) {
  using MuonId = ActsExamples::MuonSpacePoint::MuonId;

  MuonId id{};
  id.setChamber(MuonId::StationName::BIL, MuonId::DetSide::A, 1u,
                MuonId::TechField::Mdt);
  id.setLayAndCh(layer, channel);
  id.setCoordFlags(true, false, false);

  return id;
}

ActsExamples::MuonSpacePoint::MuonId makeRpcPhiId(std::uint8_t layer,
                                                  std::uint16_t channel) {
  using MuonId = ActsExamples::MuonSpacePoint::MuonId;

  MuonId id{};
  id.setChamber(MuonId::StationName::BIL, MuonId::DetSide::A, 1u,
                MuonId::TechField::Rpc);
  id.setLayAndCh(layer, channel);
  id.setCoordFlags(false, true, false);

  return id;
}

/// @brief Prepare example MuonSpaceContainer with one bucket based on test case from HoughTransformUtilsTests.cpp
ActsExamples::MuonSpacePointContainer makeDriftCircleSpacePoints(
    bool addPhiHits = false) {
  constexpr double uncert = 0.3;

  const std::array<DriftCircle, 6> driftCircles{
      DriftCircle{-427.981, -225.541, 14.5202, uncert},
      DriftCircle{-412.964, -199.530, 1.66237, uncert},
      DriftCircle{-427.981, -173.519, 12.3176, uncert},
      DriftCircle{-427.981, 173.519, 1.5412, uncert},
      DriftCircle{-442.999, 199.530, 12.3937, uncert},
      DriftCircle{-427.981, 225.541, 3.77967, uncert},
  };

  ActsExamples::MuonSpacePointContainer spacePoints{};
  spacePoints.emplace_back();

  auto& bucket = spacePoints.back();
  bucket.reserve(driftCircles.size());

  for (std::size_t i = 0; i < driftCircles.size(); ++i) {
    const DriftCircle& dc = driftCircles[i];

    ActsExamples::MuonSpacePoint& sp = bucket.emplace_back();

    sp.setGeometryId(Acts::GeometryIdentifier{i + 1u});
    sp.setId(makeMdtEtaId(static_cast<std::uint8_t>(i + 1u),
                          static_cast<std::uint16_t>(i + 1u)));

    sp.defineCoordinates(Acts::Vector3{0.0, dc.y, dc.z}, Acts::Vector3::UnitX(),
                         Acts::Vector3::UnitY());

    sp.setRadius(dc.rDrift);
    sp.setTime(0.0);
    sp.setCovariance(0.0, dc.rDriftError * dc.rDriftError, 0.0);
  }

  if (addPhiHits) {
    for (std::size_t i = 0; i < 2u; ++i) {
      ActsExamples::MuonSpacePoint& sp = bucket.emplace_back();
      sp.setGeometryId(Acts::GeometryIdentifier{driftCircles.size() + i + 1u});
      sp.setId(makeRpcPhiId(static_cast<std::uint8_t>(i + 1u),
                            static_cast<std::uint16_t>(i + 1u)));
      sp.defineCoordinates(
          Acts::Vector3{2.0 * Acts::UnitConstants::mm, 0.0,
                        static_cast<double>(i) * Acts::UnitConstants::m},
          Acts::Vector3::UnitX(), Acts::Vector3::UnitY());
      sp.setCovariance(std::pow(0.1 * Acts::UnitConstants::mm, 2), 0.0, 0.0);
    }
  }

  return spacePoints;
}

}  // namespace

BOOST_AUTO_TEST_SUITE(MuonHoughTransformSuite)

BOOST_AUTO_TEST_CASE(muon_hough_seeder_drift_circle_sanity) {
  // Truth from the original HoughTransformUtils unit test.
  constexpr double trueTanTheta = -0.0401472 / 0.994974;
  constexpr double trueInterceptY = -422.612;

  ActsExamples::MuonHoughSeeder::Config cfg{};
  cfg.inSpacePoints = "MuonSpacePoints";
  cfg.inTruthSegments =
      "TruthSegments";  // this is required even if empty in ctx
  cfg.outHoughMax = "MuonHoughMaxima";

  cfg.nBinsTanTheta = 1000;
  cfg.nBinsY0 = 1000;

  cfg.nBinsTanPhi = 10;
  cfg.nBinsX0 = 10;

  cfg.etaPlaneMarginIcept = 10.0 * Acts::UnitConstants::cm;
  cfg.phiPlaneMarginIcept = 10.0 * Acts::UnitConstants::cm;

  cfg.dumpVisualization = false;

  ActsExamples::MuonHoughSeeder seeder{
      cfg, Acts::getDefaultLogger("MuonHoughSeederTest", Acts::Logging::INFO)};

  ActsExamples::WhiteBoard eventStore{};
  ActsExamples::AlgorithmContext ctx{0, 0, eventStore, 0};

  ActsExamples::WriteDataHandle<ActsExamples::MuonSpacePointContainer>
      spacePointHandle{&seeder, "TestInputSpacePoints"};
  spacePointHandle.initialize(cfg.inSpacePoints);
  spacePointHandle(ctx, makeDriftCircleSpacePoints());

  BOOST_REQUIRE(seeder.execute(ctx) == ActsExamples::ProcessCode::SUCCESS);

  ActsExamples::ReadDataHandle<ActsExamples::MuonHoughMaxContainer>
      outputHandle{&seeder, "TestOutputHoughMaxima"};
  outputHandle.initialize(cfg.outHoughMax);

  const ActsExamples::MuonHoughMaxContainer& maxima = outputHandle(ctx);

  BOOST_REQUIRE_GT(maxima.size(), 0u);

  bool foundExpectedMaximum = false;

  double foundTanBeta = std::numeric_limits<double>::quiet_NaN();
  double foundInterceptY = std::numeric_limits<double>::quiet_NaN();
  for (const ActsExamples::MuonHoughMaximum& maximum : maxima) {
    const double tanTheta = maximum.tanBeta();
    const double interceptY = maximum.interceptY();

    if (std::abs(tanTheta - trueTanTheta) < 0.02 &&
        std::abs(interceptY - trueInterceptY) < 20.0) {
      foundExpectedMaximum = true;
      foundTanBeta = tanTheta;
      foundInterceptY = interceptY;
      break;
    }
  }

  /// Result acquired: -0.041608695652173948 -421.93644749999976
  BOOST_TEST_MESSAGE("Maximum coordinates (tanBeta, interceptY): "
                     << foundTanBeta << " " << foundInterceptY
                     << " truth: " << trueTanTheta << " " << trueInterceptY);
  BOOST_CHECK(foundExpectedMaximum);
}

BOOST_AUTO_TEST_CASE(muon_hough_seeder_drift_circle_sanity_phi_disabled) {
  // Truth from the original HoughTransformUtils unit test.
  constexpr double trueTanTheta = -0.0401472 / 0.994974;
  constexpr double trueInterceptY = -422.612;

  ActsExamples::MuonHoughSeeder::Config cfg{};
  cfg.inSpacePoints = "MuonSpacePoints";
  cfg.inTruthSegments =
      "TruthSegments";  // this is required even if empty in ctx
  cfg.outHoughMax = "MuonHoughMaxima";

  cfg.nBinsTanTheta = 1000;
  cfg.nBinsY0 = 1000;

  cfg.extendWithPhi = false;

  cfg.etaPlaneMarginIcept = 10.0 * Acts::UnitConstants::cm;
  cfg.phiPlaneMarginIcept = 10.0 * Acts::UnitConstants::cm;

  cfg.dumpVisualization = false;

  ActsExamples::MuonHoughSeeder seeder{
      cfg, Acts::getDefaultLogger("MuonHoughSeederTest", Acts::Logging::INFO)};

  ActsExamples::WhiteBoard eventStore{};
  ActsExamples::AlgorithmContext ctx{0, 0, eventStore, 0};

  ActsExamples::WriteDataHandle<ActsExamples::MuonSpacePointContainer>
      spacePointHandle{&seeder, "TestInputSpacePoints"};
  spacePointHandle.initialize(cfg.inSpacePoints);
  spacePointHandle(ctx, makeDriftCircleSpacePoints(true));

  BOOST_REQUIRE(seeder.execute(ctx) == ActsExamples::ProcessCode::SUCCESS);

  ActsExamples::ReadDataHandle<ActsExamples::MuonHoughMaxContainer>
      outputHandle{&seeder, "TestOutputHoughMaxima"};
  outputHandle.initialize(cfg.outHoughMax);

  const ActsExamples::MuonHoughMaxContainer& maxima = outputHandle(ctx);

  BOOST_REQUIRE_GT(maxima.size(), 0u);
  BOOST_CHECK(std::all_of(
      maxima.begin(), maxima.end(),
      [](const auto& maximum) { return maximum.interceptX() == 0.0; }));

  bool foundExpectedMaximum = false;

  double foundTanBeta = std::numeric_limits<double>::quiet_NaN();
  double foundInterceptY = std::numeric_limits<double>::quiet_NaN();
  for (const ActsExamples::MuonHoughMaximum& maximum : maxima) {
    const double tanTheta = maximum.tanBeta();
    const double interceptY = maximum.interceptY();

    if (std::abs(tanTheta - trueTanTheta) < 0.02 &&
        std::abs(interceptY - trueInterceptY) < 20.0) {
      foundExpectedMaximum = true;
      foundTanBeta = tanTheta;
      foundInterceptY = interceptY;
      break;
    }
  }

  /// Result acquired: -0.041608695652173948 -421.93644749999976
  BOOST_TEST_MESSAGE("Maximum coordinates (tanBeta, interceptY): "
                     << foundTanBeta << " " << foundInterceptY
                     << " truth: " << trueTanTheta << " " << trueInterceptY);
  BOOST_CHECK(foundExpectedMaximum);
}

BOOST_AUTO_TEST_CASE(muon_hough_seeder_keeps_eta_seed_when_phi_search_fails) {
  using Seeder = ActsExamples::MuonHoughSeeder;

  auto spacePoints = makeDriftCircleSpacePoints(true);
  // The two Phi hits have the same z but separated x coordinates. Their
  // Hough bands never overlap, so two input hits cannot form a Phi peak.
  spacePoints.front().back().defineCoordinates(
      Acts::Vector3{Acts::UnitConstants::m, 0.0, 0.0}, Acts::Vector3::UnitX(),
      Acts::Vector3::UnitY());

  Seeder::Config cfg{};
  cfg.inSpacePoints = "MuonSpacePoints";
  cfg.inTruthSegments = "TruthSegments";
  cfg.outHoughMax = "ExtendedMaxima";
  cfg.nBinsTanTheta = 1000;
  cfg.nBinsY0 = 1000;
  cfg.dumpVisualization = true;

  std::size_t phiSearches = 0u;
  cfg.visualizationFunction =
      [&phiSearches](
          const std::string&, const Seeder::MuonId& bucketId,
          const Seeder::MaximumVec_t& maxima, const Seeder::HoughPlane_t& plane,
          const Seeder::AxisRange_t&, const ActsExamples::MuonSegmentContainer&,
          const Acts::Logger&) {
        if (bucketId.measuresPhi()) {
          ++phiSearches;
          BOOST_CHECK_EQUAL(plane.maxHits(), 1.0f);
          BOOST_CHECK(maxima.empty());
        }
      };

  auto etaCfg = cfg;
  etaCfg.extendWithPhi = false;
  etaCfg.outHoughMax = "EtaMaxima";
  Seeder etaSeeder{
      etaCfg, Acts::getDefaultLogger("EtaSeederTest", Acts::Logging::INFO)};
  Seeder extendedSeeder{
      cfg, Acts::getDefaultLogger("ExtendedSeederTest", Acts::Logging::INFO)};

  ActsExamples::WhiteBoard eventStore{};
  ActsExamples::AlgorithmContext ctx{0u, 0u, eventStore, 0u};
  ActsExamples::WriteDataHandle<ActsExamples::MuonSpacePointContainer>
      inputHandle{&etaSeeder, "TestInputSpacePoints"};
  inputHandle.initialize(cfg.inSpacePoints);
  inputHandle(ctx, std::move(spacePoints));
  ActsExamples::WriteDataHandle<ActsExamples::MuonSegmentContainer> truthHandle{
      &etaSeeder, "TestInputTruthSegments"};
  truthHandle.initialize(cfg.inTruthSegments);
  truthHandle(ctx, ActsExamples::MuonSegmentContainer{});

  BOOST_REQUIRE(etaSeeder.execute(ctx) == ActsExamples::ProcessCode::SUCCESS);
  ActsExamples::ReadDataHandle<ActsExamples::MuonHoughMaxContainer> etaHandle{
      &etaSeeder, "TestEtaMaxima"};
  etaHandle.initialize(etaCfg.outHoughMax);
  const auto& etaMaxima = etaHandle(ctx);
  BOOST_REQUIRE(!etaMaxima.empty());

  BOOST_REQUIRE(extendedSeeder.execute(ctx) ==
                ActsExamples::ProcessCode::SUCCESS);
  ActsExamples::ReadDataHandle<ActsExamples::MuonHoughMaxContainer>
      extendedHandle{&extendedSeeder, "TestExtendedMaxima"};
  extendedHandle.initialize(cfg.outHoughMax);
  const auto& extendedMaxima = extendedHandle(ctx);

  BOOST_REQUIRE_EQUAL(phiSearches, etaMaxima.size());
  BOOST_REQUIRE_EQUAL(extendedMaxima.size(), etaMaxima.size());
  for (std::size_t i = 0u; i < etaMaxima.size(); ++i) {
    const auto& eta = etaMaxima[i];
    const auto& extended = extendedMaxima[i];
    BOOST_CHECK_EQUAL(extended.tanBeta(), eta.tanBeta());
    BOOST_CHECK_EQUAL(extended.interceptY(), eta.interceptY());
    BOOST_CHECK_EQUAL(extended.tanAlpha(), 0.0);
    BOOST_CHECK_EQUAL(extended.interceptX(), 0.0);
    for (const auto* hit : eta.hits()) {
      BOOST_CHECK(std::find(extended.hits().begin(), extended.hits().end(),
                            hit) != extended.hits().end());
    }
  }
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace ActsTests

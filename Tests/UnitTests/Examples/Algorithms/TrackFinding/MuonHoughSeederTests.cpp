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

#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
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

/// @brief Prepare example MuonSpaceContainer with one bucket based on test case from HoughTransformUtilsTests.cpp
ActsExamples::MuonSpacePointContainer makeDriftCircleSpacePoints() {
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

  cfg.etaPlaneMarginIcept = 2.0 * Acts::UnitConstants::m;
  cfg.phiPlaneMarginIcept = 2.0 * Acts::UnitConstants::m;

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

BOOST_AUTO_TEST_SUITE_END()

}  // namespace ActsTests

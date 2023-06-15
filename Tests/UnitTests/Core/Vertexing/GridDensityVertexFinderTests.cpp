// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/test/data/test_case.hpp>
#include <boost/test/tools/output_test_stream.hpp>
#include <boost/test/unit_test.hpp>

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/Common.hpp"
#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/Definitions/Units.hpp"
#include "Acts/EventData/Charge.hpp"
#include "Acts/EventData/SingleBoundTrackParameters.hpp"
#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/MagneticField/MagneticFieldContext.hpp"
#include "Acts/Surfaces/PerigeeSurface.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Tests/CommonHelpers/FloatComparisons.hpp"
#include "Acts/Utilities/Intersection.hpp"
#include "Acts/Utilities/Result.hpp"
#include "Acts/Utilities/UnitVectors.hpp"
#include "Acts/Vertexing/AdaptiveGridDensityVertexFinder.hpp"
#include "Acts/Vertexing/AdaptiveGridTrackDensity.hpp"
#include "Acts/Vertexing/GaussianGridTrackDensity.hpp"
#include "Acts/Vertexing/GridDensityVertexFinder.hpp"
#include "Acts/Vertexing/Vertex.hpp"
#include "Acts/Vertexing/VertexingOptions.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <memory>
#include <random>
#include <string>
#include <system_error>
#include <vector>

namespace bdata = boost::unit_test::data;
using namespace Acts::UnitLiterals;
using Acts::VectorHelpers::makeVector4;

namespace Acts {
namespace Test {

using Covariance = BoundSymMatrix;

// Create a test context
GeometryContext geoContext = GeometryContext();
MagneticFieldContext magFieldContext = MagneticFieldContext();

const double zVertexPos1 = 12.;
const double zVertexPos2 = -3.;
// x position
std::normal_distribution<double> xdist(1_mm, 0.1_mm);
// y position
std::normal_distribution<double> ydist(-0.7_mm, 0.1_mm);
// z1 position
std::normal_distribution<double> z1dist(zVertexPos1 * 1_mm, 1_mm);
// z2 position
std::normal_distribution<double> z2dist(zVertexPos2 * 1_mm, 0.5_mm);
// Track pT distribution
std::uniform_real_distribution<double> pTDist(0.1_GeV, 100_GeV);
// Track phi distribution
std::uniform_real_distribution<double> phiDist(-M_PI, M_PI);
// Track eta distribution
std::uniform_real_distribution<double> etaDist(-4., 4.);

///
/// @brief Unit test for GridDensityVertexFinder without caching
/// of track density values
///
BOOST_AUTO_TEST_CASE(grid_density_vertex_finder_test) {
  bool debugMode = true;

  const int mainGridSize = 3000;
  const int trkGridSize = 35;

  Covariance covMat = Covariance::Identity();

  // Perigee surface for track parameters
  Vector3 pos0{0, 0, 0};
  std::shared_ptr<PerigeeSurface> perigeeSurface =
      Surface::makeShared<PerigeeSurface>(pos0);

  VertexingOptions<BoundTrackParameters> vertexingOptions(geoContext,
                                                          magFieldContext);

  using Finder1 = GridDensityVertexFinder<mainGridSize, trkGridSize>;
  Finder1::Config cfg1;
  cfg1.cacheGridStateForTrackRemoval = false;
  Finder1 finder1(cfg1);
  Finder1::State state1;

  using AdapticeGridDensity = AdaptiveGridTrackDensity<trkGridSize>;
  // Use custom grid density here with same bin size as Finder1
  AdapticeGridDensity::Config adaptiveDensityConfig(2. / 30. * 1_mm);
  AdapticeGridDensity adaptiveDensity(adaptiveDensityConfig);

  using Finder2 = AdaptiveGridDensityVertexFinder<trkGridSize>;
  Finder2::Config cfg2(adaptiveDensity);
  cfg2.cacheGridStateForTrackRemoval = false;
  Finder2 finder2(cfg2);
  Finder2::State state2;

  int mySeed = 31415;
  std::mt19937 gen(mySeed);
  unsigned int nTracks = 200;

  std::vector<BoundTrackParameters> trackVec;
  trackVec.reserve(nTracks);

  // Create nTracks tracks for test case
  for (unsigned int i = 0; i < nTracks; i++) {
    // The position of the particle
    Vector3 pos(xdist(gen), ydist(gen), 0);

    // Create momentum and charge of track
    double pt = pTDist(gen);
    double phi = phiDist(gen);
    double eta = etaDist(gen);
    double charge = etaDist(gen) > 0 ? 1 : -1;

    // project the position on the surface
    Vector3 direction = makeDirectionUnitFromPhiEta(phi, eta);
    auto intersection = perigeeSurface->intersect(geoContext, pos, direction);
    pos = intersection.intersection.position;

    // Produce most of the tracks at near z1 position,
    // some near z2. Highest track density then expected at z1
    pos[eZ] = ((i % 4) == 0) ? z2dist(gen) : z1dist(gen);

    trackVec.push_back(BoundTrackParameters::create(
                           perigeeSurface, geoContext, makeVector4(pos, 0),
                           direction, pt, charge, covMat)
                           .value());
  }

  std::vector<const BoundTrackParameters*> trackPtrVec;
  for (const auto& trk : trackVec) {
    trackPtrVec.push_back(&trk);
  }

  auto res1 = finder1.find(trackPtrVec, vertexingOptions, state1);
  if (!res1.ok()) {
    std::cout << res1.error().message() << std::endl;
  }

  auto res2 = finder2.find(trackPtrVec, vertexingOptions, state2);
  if (!res2.ok()) {
    std::cout << res2.error().message() << std::endl;
  }

  double zResult1 = 0;
  if (res1.ok()) {
    BOOST_CHECK(!(*res1).empty());
    Vector3 result1 = (*res1).back().position();
    if (debugMode) {
      std::cout << "Vertex position result 1: " << result1 << std::endl;
    }
    CHECK_CLOSE_ABS(result1[eZ], zVertexPos1, 1_mm);
    zResult1 = result1[eZ];
  }

  double zResult2 = 0;
  if (res2.ok()) {
    BOOST_CHECK(!(*res2).empty());
    Vector3 result2 = (*res2).back().position();
    if (debugMode) {
      std::cout << "Vertex position result 2: " << result2 << std::endl;
    }
    CHECK_CLOSE_ABS(result2[eZ], zVertexPos1, 1_mm);
    zResult2 = result2[eZ];
  }

  // Both finders should give same results
  BOOST_CHECK(zResult1 == zResult2);
}

BOOST_AUTO_TEST_CASE(grid_density_vertex_finder_track_caching_test) {
  bool debugMode = true;

  const int mainGridSize = 3000;
  const int trkGridSize = 35;

  Covariance covMat = Covariance::Identity();

  // Perigee surface for track parameters
  Vector3 pos0{0, 0, 0};
  std::shared_ptr<PerigeeSurface> perigeeSurface =
      Surface::makeShared<PerigeeSurface>(pos0);

  VertexingOptions<BoundTrackParameters> vertexingOptions(geoContext,
                                                          magFieldContext);

  using Finder1 = GridDensityVertexFinder<mainGridSize, trkGridSize>;
  using GridDensity = GaussianGridTrackDensity<mainGridSize, trkGridSize>;

  // Use custom grid density here
  GridDensity::Config densityConfig(100_mm);
  densityConfig.useHighestSumZPosition = true;
  GridDensity density(densityConfig);

  Finder1::Config cfg(density);
  cfg.cacheGridStateForTrackRemoval = true;
  Finder1 finder1(cfg);

  using AdapticeGridDensity = AdaptiveGridTrackDensity<trkGridSize>;
  // Use custom grid density here with same bin size as Finder1
  AdapticeGridDensity::Config adaptiveDensityConfig(2. / 30. * 1_mm);
  adaptiveDensityConfig.useHighestSumZPosition = true;
  AdapticeGridDensity adaptiveDensity(adaptiveDensityConfig);

  using Finder2 = AdaptiveGridDensityVertexFinder<trkGridSize>;
  Finder2::Config cfg2(adaptiveDensity);
  cfg2.cacheGridStateForTrackRemoval = true;
  Finder2 finder2(cfg2);

  int mySeed = 31415;
  std::mt19937 gen(mySeed);
  unsigned int nTracks = 200;

  std::vector<BoundTrackParameters> trackVec;
  trackVec.reserve(nTracks);

  // Create nTracks tracks for test case
  for (unsigned int i = 0; i < nTracks; i++) {
    // The position of the particle
    Vector3 pos(xdist(gen), ydist(gen), 0);

    // Create momentum and charge of track
    double pt = pTDist(gen);
    double phi = phiDist(gen);
    double eta = etaDist(gen);
    double charge = etaDist(gen) > 0 ? 1 : -1;

    // project the position on the surface
    Vector3 direction = makeDirectionUnitFromPhiEta(phi, eta);
    auto intersection = perigeeSurface->intersect(geoContext, pos, direction);
    pos = intersection.intersection.position;

    // Produce most of the tracks at near z1 position,
    // some near z2. Highest track density then expected at z1
    pos[eZ] = ((i % 4) == 0) ? z2dist(gen) : z1dist(gen);

    trackVec.push_back(BoundTrackParameters::create(
                           perigeeSurface, geoContext, makeVector4(pos, 0),
                           direction, pt, charge, covMat)
                           .value());
  }

  std::vector<const BoundTrackParameters*> trackPtrVec;
  for (const auto& trk : trackVec) {
    trackPtrVec.push_back(&trk);
  }

  Finder1::State state1;
  Finder2::State state2;

  double zResult1 = 0;
  double zResult2 = 0;

  auto res1 = finder1.find(trackPtrVec, vertexingOptions, state1);
  if (!res1.ok()) {
    std::cout << res1.error().message() << std::endl;
  }
  if (res1.ok()) {
    BOOST_CHECK(!(*res1).empty());
    Vector3 result = (*res1).back().position();
    if (debugMode) {
      std::cout << "Vertex position after first fill 1: " << result
                << std::endl;
    }
    CHECK_CLOSE_ABS(result[eZ], zVertexPos1, 1_mm);
    zResult1 = result[eZ];
  }

  auto res2 = finder2.find(trackPtrVec, vertexingOptions, state2);
  if (!res2.ok()) {
    std::cout << res2.error().message() << std::endl;
  }
  if (res2.ok()) {
    BOOST_CHECK(!(*res2).empty());
    Vector3 result = (*res2).back().position();
    if (debugMode) {
      std::cout << "Vertex position after first fill 2: " << result
                << std::endl;
    }
    CHECK_CLOSE_ABS(result[eZ], zVertexPos1, 1_mm);
    zResult2 = result[eZ];
  }

  BOOST_CHECK(zResult1 == zResult2);

  int trkCount = 0;
  std::vector<const BoundTrackParameters*> removedTracks;
  for (const auto& trk : trackVec) {
    if ((trkCount % 4) != 0) {
      removedTracks.push_back(&trk);
    }
    trkCount++;
  }

  state1.tracksToRemove = removedTracks;
  state2.tracksToRemove = removedTracks;

  auto res3 = finder1.find(trackPtrVec, vertexingOptions, state1);
  if (!res3.ok()) {
    std::cout << res3.error().message() << std::endl;
  }
  if (res3.ok()) {
    BOOST_CHECK(!(*res3).empty());
    Vector3 result = (*res3).back().position();
    if (debugMode) {
      std::cout
          << "Vertex position after removing tracks near first density peak 1: "
          << result << std::endl;
    }
    CHECK_CLOSE_ABS(result[eZ], zVertexPos2, 1_mm);
    zResult1 = result[eZ];
  }

  auto res4 = finder2.find(trackPtrVec, vertexingOptions, state2);
  if (!res4.ok()) {
    std::cout << res4.error().message() << std::endl;
  }
  if (res4.ok()) {
    BOOST_CHECK(!(*res4).empty());
    Vector3 result = (*res4).back().position();
    if (debugMode) {
      std::cout
          << "Vertex position after removing tracks near first density peak 2: "
          << result << std::endl;
    }
    CHECK_CLOSE_ABS(result[eZ], zVertexPos2, 1_mm);
    zResult2 = result[eZ];
  }

  BOOST_CHECK(zResult1 == zResult2);
}

///
/// @brief Unit test for GridDensityVertexFinder with seed with estimation
///
BOOST_AUTO_TEST_CASE(grid_density_vertex_finder_seed_width_test) {
  bool debugMode = true;

  const int mainGridSize = 3000;
  const int trkGridSize = 35;

  Covariance covMat = Covariance::Identity();

  // Perigee surface for track parameters
  Vector3 pos0{0, 0, 0};
  std::shared_ptr<PerigeeSurface> perigeeSurface =
      Surface::makeShared<PerigeeSurface>(pos0);

  VertexingOptions<BoundTrackParameters> vertexingOptions(geoContext,
                                                          magFieldContext);
  Vertex<BoundTrackParameters> constraintVtx;
  constraintVtx.setCovariance(SymMatrix3::Identity());
  vertexingOptions.vertexConstraint = constraintVtx;

  using Finder1 = GridDensityVertexFinder<mainGridSize, trkGridSize>;
  Finder1::Config cfg1;
  cfg1.cacheGridStateForTrackRemoval = false;
  cfg1.estimateSeedWidth = true;
  Finder1 finder1(cfg1);
  Finder1::State state1;

  using AdapticeGridDensity = AdaptiveGridTrackDensity<trkGridSize>;
  // Use custom grid density here with same bin size as Finder1
  AdapticeGridDensity::Config adaptiveDensityConfig(2. / 30. * 1_mm);
  AdapticeGridDensity adaptiveDensity(adaptiveDensityConfig);

  using Finder2 = AdaptiveGridDensityVertexFinder<trkGridSize>;
  Finder2::Config cfg2(adaptiveDensity);
  cfg2.cacheGridStateForTrackRemoval = false;
  cfg2.estimateSeedWidth = true;
  Finder2 finder2(cfg2);
  Finder2::State state2;

  int mySeed = 31415;
  std::mt19937 gen(mySeed);
  unsigned int nTracks = 20;

  std::vector<BoundTrackParameters> trackVec;
  trackVec.reserve(nTracks);

  // Create nTracks tracks for test case
  for (unsigned int i = 0; i < nTracks; i++) {
    // The position of the particle
    Vector3 pos(xdist(gen), ydist(gen), 0);

    // Create momentum and charge of track
    double pt = pTDist(gen);
    double phi = phiDist(gen);
    double eta = etaDist(gen);
    double charge = etaDist(gen) > 0 ? 1 : -1;

    // project the position on the surface
    Vector3 direction = makeDirectionUnitFromPhiEta(phi, eta);
    auto intersection = perigeeSurface->intersect(geoContext, pos, direction);
    pos = intersection.intersection.position;

    pos[eZ] = z1dist(gen);

    trackVec.push_back(BoundTrackParameters::create(
                           perigeeSurface, geoContext, makeVector4(pos, 0),
                           direction, pt, charge, covMat)
                           .value());
  }

  std::vector<const BoundTrackParameters*> trackPtrVec;
  for (const auto& trk : trackVec) {
    trackPtrVec.push_back(&trk);
  }

  // Test finder 1
  auto res1 = finder1.find(trackPtrVec, vertexingOptions, state1);
  if (!res1.ok()) {
    std::cout << res1.error().message() << std::endl;
  }

  double covZZ1 = 0;
  if (res1.ok()) {
    BOOST_CHECK(!(*res1).empty());
    SymMatrix3 cov = (*res1).back().covariance();
    BOOST_CHECK(constraintVtx.covariance() != cov);
    BOOST_CHECK(cov(eZ, eZ) != 0.);
    covZZ1 = cov(eZ, eZ);
    if (debugMode) {
      std::cout << "Estimated z-seed width 1: " << cov(eZ, eZ) << std::endl;
    }
  }

  // Test finder 2
  auto res2 = finder2.find(trackPtrVec, vertexingOptions, state2);
  if (!res2.ok()) {
    std::cout << res2.error().message() << std::endl;
  }

  double covZZ2 = 0;
  if (res2.ok()) {
    BOOST_CHECK(!(*res2).empty());
    SymMatrix3 cov = (*res2).back().covariance();
    BOOST_CHECK(constraintVtx.covariance() != cov);
    BOOST_CHECK(cov(eZ, eZ) != 0.);
    covZZ2 = cov(eZ, eZ);
    if (debugMode) {
      std::cout << "Estimated z-seed width 2: " << cov(eZ, eZ) << std::endl;
    }
  }

  // Test for same seed width
  CHECK_CLOSE_ABS(covZZ1, covZZ2, 1e-4);
}

}  // namespace Test
}  // namespace Acts

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

#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/Surfaces/PerigeeSurface.hpp"
#include "Acts/Tests/CommonHelpers/FloatComparisons.hpp"
#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Utilities/Units.hpp"
#include "Acts/Vertexing/GaussianGridTrackDensity.hpp"
#include "Acts/Vertexing/GridDensityVertexFinder.hpp"

namespace bdata = boost::unit_test::data;
using namespace Acts::UnitLiterals;

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
  Vector3D pos0{0, 0, 0};
  std::shared_ptr<PerigeeSurface> perigeeSurface =
      Surface::makeShared<PerigeeSurface>(pos0);

  VertexingOptions<BoundParameters> vertexingOptions(geoContext,
                                                     magFieldContext);

  using Finder = GridDensityVertexFinder<mainGridSize, trkGridSize>;
  Finder::Config cfg;
  cfg.cacheGridStateForTrackRemoval = false;
  Finder finder(cfg);

  Finder::State state;

  int mySeed = 31415;
  std::mt19937 gen(mySeed);
  unsigned int nTracks = 200;

  std::vector<BoundParameters> trackVec;
  trackVec.reserve(nTracks);

  // Create nTracks tracks for test case
  for (unsigned int i = 0; i < nTracks; i++) {
    // The position of the particle
    Vector3D pos(xdist(gen), ydist(gen), 0);
    // Produce most of the tracks at near z1 position,
    // some near z2. Highest track density then expected at z1
    if ((i % 4) == 0) {
      pos[eZ] = z2dist(gen);
    } else {
      pos[eZ] = z1dist(gen);
    }

    // Create momentum and charge of track
    double pt = pTDist(gen);
    double phi = phiDist(gen);
    double eta = etaDist(gen);
    Vector3D mom(pt * std::cos(phi), pt * std::sin(phi), pt * std::sinh(eta));
    double charge = etaDist(gen) > 0 ? 1 : -1;

    trackVec.push_back(BoundParameters(geoContext, covMat, pos, mom, charge, 0,
                                       perigeeSurface));
  }

  std::vector<const BoundParameters*> trackPtrVec;
  for (const auto& trk : trackVec) {
    trackPtrVec.push_back(&trk);
  }

  auto res = finder.find(trackPtrVec, vertexingOptions, state);
  if (!res.ok()) {
    std::cout << res.error().message() << std::endl;
  }

  if (res.ok()) {
    BOOST_CHECK(!(*res).empty());
    Vector3D result = (*res).back().position();
    if (debugMode) {
      std::cout << "Vertex position result: " << result << std::endl;
    }
    CHECK_CLOSE_ABS(result[eZ], zVertexPos1, 1_mm);
  }
}

BOOST_AUTO_TEST_CASE(grid_density_vertex_finder_track_caching_test) {
  bool debugMode = true;

  const int mainGridSize = 3000;
  const int trkGridSize = 35;

  Covariance covMat = Covariance::Identity();

  // Perigee surface for track parameters
  Vector3D pos0{0, 0, 0};
  std::shared_ptr<PerigeeSurface> perigeeSurface =
      Surface::makeShared<PerigeeSurface>(pos0);

  VertexingOptions<BoundParameters> vertexingOptions(geoContext,
                                                     magFieldContext);

  using Finder = GridDensityVertexFinder<mainGridSize, trkGridSize>;
  using GridDensity = GaussianGridTrackDensity<mainGridSize, trkGridSize>;

  // Use custom grid density here
  GridDensity::Config densityConfig(100_mm);
  densityConfig.useHighestSumZPosition = true;
  GridDensity density(densityConfig);

  Finder::Config cfg(density);
  cfg.cacheGridStateForTrackRemoval = true;
  Finder finder(cfg);

  int mySeed = 31415;
  std::mt19937 gen(mySeed);
  unsigned int nTracks = 200;

  std::vector<BoundParameters> trackVec;
  trackVec.reserve(nTracks);

  // Create nTracks tracks for test case
  for (unsigned int i = 0; i < nTracks; i++) {
    // The position of the particle
    Vector3D pos(xdist(gen), ydist(gen), 0);
    // Produce most of the tracks at near z1 position,
    // some near z2. Highest track density then expected at z1
    if ((i % 4) == 0) {
      pos[eZ] = z2dist(gen);
    } else {
      pos[eZ] = z1dist(gen);
    }

    // Create momentum and charge of track
    double pt = pTDist(gen);
    double phi = phiDist(gen);
    double eta = etaDist(gen);
    Vector3D mom(pt * std::cos(phi), pt * std::sin(phi), pt * std::sinh(eta));
    double charge = etaDist(gen) > 0 ? 1 : -1;

    trackVec.push_back(BoundParameters(geoContext, covMat, pos, mom, charge, 0,
                                       perigeeSurface));
  }

  std::vector<const BoundParameters*> trackPtrVec;
  for (const auto& trk : trackVec) {
    trackPtrVec.push_back(&trk);
  }

  Finder::State state;

  auto res = finder.find(trackPtrVec, vertexingOptions, state);
  if (!res.ok()) {
    std::cout << res.error().message() << std::endl;
  }
  if (res.ok()) {
    BOOST_CHECK(!(*res).empty());
    Vector3D result = (*res).back().position();
    if (debugMode) {
      std::cout << "Vertex position after first fill: " << result << std::endl;
    }
    CHECK_CLOSE_ABS(result[eZ], zVertexPos1, 1_mm);
  }

  int trkCount = 0;
  std::vector<const BoundParameters*> removedTracks;
  for (const auto& trk : trackVec) {
    if ((trkCount % 4) != 0) {
      removedTracks.push_back(&trk);
    }
    trkCount++;
  }

  state.tracksToRemove = removedTracks;

  auto res2 = finder.find(trackPtrVec, vertexingOptions, state);
  if (!res2.ok()) {
    std::cout << res2.error().message() << std::endl;
  }
  if (res2.ok()) {
    BOOST_CHECK(!(*res2).empty());
    Vector3D result = (*res2).back().position();
    if (debugMode) {
      std::cout
          << "Vertex position after removing tracks near first density peak: "
          << result << std::endl;
    }
    CHECK_CLOSE_ABS(result[eZ], zVertexPos2, 1_mm);
  }
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
  Vector3D pos0{0, 0, 0};
  std::shared_ptr<PerigeeSurface> perigeeSurface =
      Surface::makeShared<PerigeeSurface>(pos0);

  VertexingOptions<BoundParameters> vertexingOptions(geoContext,
                                                     magFieldContext);
  Vertex<BoundParameters> constraintVtx;
  constraintVtx.setCovariance(ActsSymMatrixD<3>::Identity());
  vertexingOptions.vertexConstraint = constraintVtx;

  using Finder = GridDensityVertexFinder<mainGridSize, trkGridSize>;
  Finder::Config cfg;
  cfg.cacheGridStateForTrackRemoval = false;
  cfg.estimateSeedWidth = true;
  Finder finder(cfg);

  Finder::State state;

  int mySeed = 31415;
  std::mt19937 gen(mySeed);
  unsigned int nTracks = 20;

  std::vector<BoundParameters> trackVec;
  trackVec.reserve(nTracks);

  // Create nTracks tracks for test case
  for (unsigned int i = 0; i < nTracks; i++) {
    // The position of the particle
    Vector3D pos(xdist(gen), ydist(gen), z1dist(gen));

    // Create momentum and charge of track
    double pt = pTDist(gen);
    double phi = phiDist(gen);
    double eta = etaDist(gen);
    Vector3D mom(pt * std::cos(phi), pt * std::sin(phi), pt * std::sinh(eta));
    double charge = etaDist(gen) > 0 ? 1 : -1;

    trackVec.push_back(BoundParameters(geoContext, covMat, pos, mom, charge, 0,
                                       perigeeSurface));
  }

  std::vector<const BoundParameters*> trackPtrVec;
  for (const auto& trk : trackVec) {
    trackPtrVec.push_back(&trk);
  }

  auto res = finder.find(trackPtrVec, vertexingOptions, state);
  if (!res.ok()) {
    std::cout << res.error().message() << std::endl;
  }

  if (res.ok()) {
    BOOST_CHECK(!(*res).empty());
    ActsSymMatrixD<3> cov = (*res).back().covariance();
    BOOST_CHECK(constraintVtx.covariance() != cov);
    BOOST_CHECK(cov(eZ, eZ) != 0.);
    if (debugMode) {
      std::cout << "Estimated z-seed width: " << cov(eZ, eZ) << std::endl;
    }
  }
}

}  // namespace Test
}  // namespace Acts

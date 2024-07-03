// This file is part of the Acts project.
//
// Copyright (C) 2020-2024 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/test/data/test_case.hpp>
#include <boost/test/tools/output_test_stream.hpp>
#include <boost/test/unit_test.hpp>

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/Definitions/Units.hpp"
#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/Surfaces/PerigeeSurface.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Tests/CommonHelpers/FloatComparisons.hpp"
#include "Acts/Utilities/Result.hpp"
#include "Acts/Vertexing/AdaptiveGridTrackDensity.hpp"

#include <algorithm>
#include <iterator>
#include <memory>
#include <optional>
#include <utility>
#include <vector>

namespace bdata = boost::unit_test::data;
using namespace Acts::UnitLiterals;

namespace Acts::Test {

BOOST_AUTO_TEST_CASE(adaptive_gaussian_grid_density_track_adding_test) {
  const int trkGridSize = 15;

  double binSize = 0.1;  // mm

  // Set up grid density with zMinMax
  AdaptiveGridTrackDensity<trkGridSize>::Config cfg;
  cfg.binSize = binSize;
  AdaptiveGridTrackDensity<trkGridSize> grid(cfg);

  // Create some test tracks in such a way that some tracks
  //  e.g. overlap and that certain tracks need to be inserted
  // between two other tracks
  BoundSquareMatrix covMat = BoundSquareMatrix::Identity();

  BoundVector paramVec0;
  paramVec0 << 100.0, -0.45, 0, 0, 0, 0;
  BoundVector paramVec1;
  paramVec1 << 0.01, -0.45, 0, 0, 0, 0;
  BoundVector paramVec2;
  paramVec2 << 0.01, 10.95, 0, 0, 0, 0;
  BoundVector paramVec3;
  paramVec3 << 0.01, 0.95, 0, 0, 0, 0;
  BoundVector paramVec4;
  paramVec4 << 0.01, -30.95, 0, 0, 0, 0;
  BoundVector paramVec5;
  paramVec5 << 0.01, -15.0, 0, 0, 0, 0;

  // Create perigee surface
  std::shared_ptr<PerigeeSurface> perigeeSurface =
      Surface::makeShared<PerigeeSurface>(Vector3(0., 0., 0.));

  BoundTrackParameters params0(perigeeSurface, paramVec0, covMat,
                               ParticleHypothesis::pion());
  BoundTrackParameters params1(perigeeSurface, paramVec1, covMat,
                               ParticleHypothesis::pion());
  BoundTrackParameters params2(perigeeSurface, paramVec2, covMat,
                               ParticleHypothesis::pion());
  BoundTrackParameters params3(perigeeSurface, paramVec3, covMat,
                               ParticleHypothesis::pion());
  BoundTrackParameters params4(perigeeSurface, paramVec4, covMat,
                               ParticleHypothesis::pion());
  BoundTrackParameters params5(perigeeSurface, paramVec5, covMat,
                               ParticleHypothesis::pion());

  // Start with empty grids
  AdaptiveGridTrackDensity<trkGridSize>::MainDensityMap mainDensityMap;

  // Track is too far away from z axis and was not added
  auto zBinAndTrack = grid.addTrack(params0, mainDensityMap);
  BOOST_CHECK(mainDensityMap.empty());
  BOOST_CHECK_EQUAL(mainDensityMap.density.size(), mainDensityMap.zBin.size());

  // Track should have been entirely added to both grids
  zBinAndTrack = grid.addTrack(params1, mainDensityMap);
  BOOST_CHECK_EQUAL(mainDensityMap.size(), trkGridSize);
  BOOST_CHECK_EQUAL(mainDensityMap.density.size(), mainDensityMap.zBin.size());

  // Track should have been entirely added to both grids
  zBinAndTrack = grid.addTrack(params2, mainDensityMap);
  BOOST_CHECK_EQUAL(mainDensityMap.density.size(), 2 * trkGridSize);
  BOOST_CHECK_EQUAL(mainDensityMap.density.size(), mainDensityMap.zBin.size());

  // Track 3 has overlap of 2 bins with track 1
  zBinAndTrack = grid.addTrack(params3, mainDensityMap);
  BOOST_CHECK_EQUAL(mainDensityMap.density.size(), 3 * trkGridSize - 2);
  BOOST_CHECK_EQUAL(mainDensityMap.density.size(), mainDensityMap.zBin.size());

  // Add first track again, should *not* introduce new z entries
  zBinAndTrack = grid.addTrack(params1, mainDensityMap);
  BOOST_CHECK_EQUAL(mainDensityMap.density.size(), 3 * trkGridSize - 2);
  BOOST_CHECK_EQUAL(mainDensityMap.density.size(), mainDensityMap.zBin.size());

  // Add two more tracks and check if order is correct
  zBinAndTrack = grid.addTrack(params4, mainDensityMap);
  BOOST_CHECK_EQUAL(mainDensityMap.density.size(), mainDensityMap.zBin.size());
  zBinAndTrack = grid.addTrack(params5, mainDensityMap);
  BOOST_CHECK_EQUAL(mainDensityMap.density.size(), mainDensityMap.zBin.size());

  BOOST_CHECK(std::is_sorted(std::begin(mainDensityMap.zBin),
                             std::end(mainDensityMap.zBin)));
}

BOOST_AUTO_TEST_CASE(adaptive_gaussian_grid_density_max_z_and_width_test) {
  const int trkGridSize = 15;

  double binSize = 0.1;  // mm

  // Set up grid density with zMinMax
  AdaptiveGridTrackDensity<trkGridSize>::Config cfg;
  cfg.binSize = binSize;
  AdaptiveGridTrackDensity<trkGridSize> grid(cfg);

  // Create some test tracks
  BoundSquareMatrix covMat = BoundSquareMatrix::Identity();

  float z0Trk1 = 0.25;
  float z0Trk2 = -10.95;
  BoundVector paramVec1;
  paramVec1 << 0.02, z0Trk1, 0, 0, 0, 0;
  BoundVector paramVec2;
  paramVec2 << 0.01, z0Trk2, 0, 0, 0, 0;

  // Create perigee surface
  std::shared_ptr<PerigeeSurface> perigeeSurface =
      Surface::makeShared<PerigeeSurface>(Vector3(0., 0., 0.));

  BoundTrackParameters params1(perigeeSurface, paramVec1, covMat,
                               ParticleHypothesis::pion());
  BoundTrackParameters params2(perigeeSurface, paramVec2, covMat,
                               ParticleHypothesis::pion());

  // Start with empty grids
  AdaptiveGridTrackDensity<trkGridSize>::MainDensityMap mainDensityMap;

  // Fill grid with track densities
  auto zBinAndTrack = grid.addTrack(params1, mainDensityMap);
  auto res1 = grid.getMaxZPosition(mainDensityMap);
  BOOST_CHECK(res1.ok());
  // Maximum should be at z0Trk1 position
  BOOST_CHECK_EQUAL(*res1, z0Trk1);

  // Add second track
  zBinAndTrack = grid.addTrack(params2, mainDensityMap);
  auto res2 = grid.getMaxZPosition(mainDensityMap);
  BOOST_CHECK(res2.ok());
  // Trk 2 is closer to z-axis and should yield higher density values
  // New maximum is therefore at z0Trk2
  BOOST_CHECK_EQUAL(*res2, z0Trk2);

  // Get max position and width estimation
  auto resWidth1 = grid.getMaxZPositionAndWidth(mainDensityMap);
  BOOST_CHECK(resWidth1.ok());
  BOOST_CHECK_EQUAL((*resWidth1).first, z0Trk2);
  BOOST_CHECK((*resWidth1).second > 0);
}

BOOST_AUTO_TEST_CASE(adaptive_gaussian_grid_density_highest_density_sum_test) {
  const int trkGridSize = 15;

  double binSize = 0.1;  // mm

  // Set up grid density with zMinMax
  AdaptiveGridTrackDensity<trkGridSize>::Config cfg;
  cfg.binSize = binSize;
  cfg.useHighestSumZPosition = true;

  AdaptiveGridTrackDensity<trkGridSize> grid(cfg);

  // Create some test tracks
  BoundSquareMatrix covMat = BoundSquareMatrix::Identity();

  float z0Trk1 = 0.25;
  float z0Trk2 = -10.95;
  BoundVector paramVec1;
  paramVec1 << 0.01, z0Trk1, 0, 0, 0, 0;
  BoundVector paramVec2;
  paramVec2 << 0.009, z0Trk2, 0, 0, 0, 0;

  // Create perigee surface
  std::shared_ptr<PerigeeSurface> perigeeSurface =
      Surface::makeShared<PerigeeSurface>(Vector3(0., 0., 0.));

  BoundTrackParameters params1(perigeeSurface, paramVec1, covMat,
                               ParticleHypothesis::pion());
  BoundTrackParameters params2(perigeeSurface, paramVec2, covMat,
                               ParticleHypothesis::pion());

  // Start with empty grids
  AdaptiveGridTrackDensity<trkGridSize>::MainDensityMap mainDensityMap;

  // Fill grid with track densities
  auto zBinAndTrack = grid.addTrack(params1, mainDensityMap);

  auto res1 = grid.getMaxZPosition(mainDensityMap);
  BOOST_CHECK(res1.ok());
  // Maximum should be at z0Trk1 position
  BOOST_CHECK_EQUAL(*res1, z0Trk1);

  // Add second track
  zBinAndTrack = grid.addTrack(params2, mainDensityMap);
  auto res2 = grid.getMaxZPosition(mainDensityMap);
  BOOST_CHECK(res2.ok());
  // Trk 2 is closer to z-axis and should yield higher density values
  // New maximum is therefore at z0Trk2
  BOOST_CHECK_EQUAL(*res2, z0Trk2);

  // Add small density values around the maximum of track 1
  const float densityToAdd = 5e-4;
  mainDensityMap.density[21] += densityToAdd;
  mainDensityMap.density[23] += densityToAdd;

  auto res3 = grid.getMaxZPosition(mainDensityMap);
  BOOST_CHECK(res3.ok());
  // Trk 2 still has the highest peak density value, however, the small
  // added densities for track 1 around its maximum should now lead to
  // a prediction at z0Trk1 again
  BOOST_CHECK_EQUAL(*res3, z0Trk1);
}

BOOST_AUTO_TEST_CASE(adaptive_gaussian_grid_density_track_removing_test) {
  const int trkGridSize = 15;

  double binSize = 0.1;  // mm

  // Set up grid density with zMinMax
  AdaptiveGridTrackDensity<trkGridSize>::Config cfg;
  cfg.binSize = binSize;
  AdaptiveGridTrackDensity<trkGridSize> grid(cfg);

  // Create some test tracks
  BoundSquareMatrix covMat = BoundSquareMatrix::Identity();

  // Define z0 values for test tracks
  float z0Trk1 = -0.45;
  float z0Trk2 = -0.25;

  BoundVector paramVec0;
  paramVec0 << 0.1, z0Trk1, 0, 0, 0, 0;
  BoundVector paramVec1;
  paramVec1 << 0.1, z0Trk2, 0, 0, 0, 0;

  // Create perigee surface
  std::shared_ptr<PerigeeSurface> perigeeSurface =
      Surface::makeShared<PerigeeSurface>(Vector3(0., 0., 0.));

  BoundTrackParameters params0(perigeeSurface, paramVec0, covMat,
                               ParticleHypothesis::pion());
  BoundTrackParameters params1(perigeeSurface, paramVec1, covMat,
                               ParticleHypothesis::pion());

  // Start with empty grids
  AdaptiveGridTrackDensity<trkGridSize>::MainDensityMap mainDensityMap;

  // Add track 0
  auto zBinAndTrack0 = grid.addTrack(params0, mainDensityMap);
  BOOST_CHECK(!mainDensityMap.empty());
  BOOST_CHECK_EQUAL(mainDensityMap.density.size(), mainDensityMap.zBin.size());
  // Grid size should match trkGridSize
  BOOST_CHECK_EQUAL(mainDensityMap.size(), trkGridSize);

  // Calculate total density
  float densitySum0 = 0;
  for (double d : mainDensityMap.density) {
    densitySum0 += d;
  }

  // Add track 0 again
  auto zBinAndTrack1 = grid.addTrack(params0, mainDensityMap);
  BOOST_CHECK(!mainDensityMap.empty());
  BOOST_CHECK_EQUAL(mainDensityMap.density.size(), mainDensityMap.zBin.size());
  // Grid size should still match trkGridSize
  BOOST_CHECK_EQUAL(mainDensityMap.size(), trkGridSize);

  // Calculate new total density
  float densitySum1 = 0;
  for (double d : mainDensityMap.density) {
    densitySum1 += d;
  }

  BOOST_CHECK(2 * densitySum0 == densitySum1);

  // Remove track 1
  grid.subtractTrack(zBinAndTrack1, mainDensityMap);

  // Calculate new total density
  float densitySum2 = 0;
  for (double d : mainDensityMap.density) {
    densitySum2 += d;
  }

  // Density should be old one again
  BOOST_CHECK(densitySum0 == densitySum2);
  // Grid size should still match trkGridSize (removal does not touch grid size)
  BOOST_CHECK_EQUAL(mainDensityMap.size(), trkGridSize);

  // Add track 1, overlapping track 0
  auto zBinAndTrack2 = grid.addTrack(params1, mainDensityMap);
  BOOST_CHECK_EQUAL(mainDensityMap.density.size(), mainDensityMap.zBin.size());

  int nNonOverlappingBins =
      static_cast<int>(std::abs(z0Trk1 - z0Trk2) / binSize + 1);
  BOOST_CHECK_EQUAL(mainDensityMap.size(), trkGridSize + nNonOverlappingBins);

  float densitySum3 = 0;
  for (double d : mainDensityMap.density) {
    densitySum3 += d;
  }

  // Remove second track 1
  grid.subtractTrack(zBinAndTrack0, mainDensityMap);

  float densitySum4 = 0;
  for (double d : mainDensityMap.density) {
    densitySum4 += d;
  }

  // Density should match differences of removed tracks
  CHECK_CLOSE_ABS(densitySum4, densitySum3 - densitySum0, 1e-5);

  // Remove last track again
  grid.subtractTrack(zBinAndTrack2, mainDensityMap);

  // Size should not have changed
  BOOST_CHECK_EQUAL(mainDensityMap.size(), trkGridSize + nNonOverlappingBins);

  float densitySum5 = 0;
  for (double d : mainDensityMap.density) {
    densitySum5 += d;
  }

  // Grid is now empty after all tracks were removed
  CHECK_CLOSE_ABS(densitySum5, 0., 1e-5);
}

}  // namespace Acts::Test

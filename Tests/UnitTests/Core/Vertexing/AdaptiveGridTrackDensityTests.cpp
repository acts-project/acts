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
#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/Definitions/Units.hpp"
#include "Acts/EventData/GenericBoundTrackParameters.hpp"
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

namespace Acts {
namespace Test {

using Covariance = BoundSymMatrix;

BOOST_AUTO_TEST_CASE(adaptive_gaussian_grid_density_track_adding_test) {
  const int trkGridSize = 15;

  double binSize = 0.1;  // mm

  // Set up grid density with zMinMax
  AdaptiveGridTrackDensity<trkGridSize>::Config cfg(binSize);
  AdaptiveGridTrackDensity<trkGridSize> grid(cfg);

  // Create some test tracks in such a way that some tracks
  //  e.g. overlap and that certain tracks need to be inserted
  // between two other tracks
  Covariance covMat(Covariance::Identity());

  BoundVector paramVec0;
  paramVec0 << 100.0, -0.45, 0, 0, 0, 0;
  BoundVector paramVec1;
  paramVec1 << 0.01, -0.45, 0, 0, 0, 0;
  BoundVector paramVec2;
  paramVec2 << 0.01, 10.95, 0, 0, 0, 0;
  BoundVector paramVec3;
  paramVec3 << 0.01, 0.95, 0, 0, 0, 0;

  // Create perigee surface
  std::shared_ptr<PerigeeSurface> perigeeSurface =
      Surface::makeShared<PerigeeSurface>(Vector3(0., 0., 0.));

  BoundTrackParameters params0(perigeeSurface, paramVec0, covMat);
  BoundTrackParameters params1(perigeeSurface, paramVec1, covMat);
  BoundTrackParameters params2(perigeeSurface, paramVec2, covMat);
  BoundTrackParameters params3(perigeeSurface, paramVec3, covMat);

  // Empty map
  AdaptiveGridTrackDensity<trkGridSize>::DensityMap mainDensityMap;

  // Track is too far away from z axis and was not added
  auto trackDensityMap = grid.addTrack(params0, mainDensityMap);
  BOOST_CHECK(mainDensityMap.empty());

  // Track should have been entirely added to both grids
  trackDensityMap = grid.addTrack(params1, mainDensityMap);
  BOOST_CHECK_EQUAL(mainDensityMap.size(), trkGridSize);

  // Track should have been entirely added to both grids
  trackDensityMap = grid.addTrack(params2, mainDensityMap);
  BOOST_CHECK_EQUAL(mainDensityMap.size(), 2 * trkGridSize);

  // Track 3 has overlap of 2 bins with track 1
  trackDensityMap = grid.addTrack(params3, mainDensityMap);
  BOOST_CHECK_EQUAL(mainDensityMap.size(), 3 * trkGridSize - 2);

  // Add first track again, should *not* introduce new z entries
  trackDensityMap = grid.addTrack(params1, mainDensityMap);
  BOOST_CHECK_EQUAL(mainDensityMap.size(), 3 * trkGridSize - 2);
}

BOOST_AUTO_TEST_CASE(adaptive_gaussian_grid_density_max_z_and_width_test) {
  const int trkGridSize = 15;

  double binSize = 0.1;  // mm

  // Set up grid density with zMinMax
  AdaptiveGridTrackDensity<trkGridSize>::Config cfg(binSize);
  AdaptiveGridTrackDensity<trkGridSize> grid(cfg);

  // Create some test tracks
  Covariance covMat(Covariance::Identity() * 0.005);

  float z0Trk1 = 0.25;
  float z0Trk2 = -10.95;
  BoundVector paramVec1;
  paramVec1 << 0.02, z0Trk1, 0, 0, 0, 0;
  BoundVector paramVec2;
  paramVec2 << 0.01, z0Trk2, 0, 0, 0, 0;

  // Create perigee surface
  std::shared_ptr<PerigeeSurface> perigeeSurface =
      Surface::makeShared<PerigeeSurface>(Vector3(0., 0., 0.));

  BoundTrackParameters params1(perigeeSurface, paramVec1, covMat);
  BoundTrackParameters params2(perigeeSurface, paramVec2, covMat);

  // Empty map
  AdaptiveGridTrackDensity<trkGridSize>::DensityMap mainDensityMap;

  auto trackDensityMap = grid.addTrack(params1, mainDensityMap);
  auto res1 = grid.getMaxZPosition(mainDensityMap);
  BOOST_CHECK(res1.ok());
  // Maximum should be at z0Trk1 position
  BOOST_CHECK_EQUAL(*res1, z0Trk1);

  trackDensityMap = grid.addTrack(params2, mainDensityMap);
  auto res2 = grid.getMaxZPosition(mainDensityMap);
  BOOST_CHECK(res2.ok());
  // Trk 2 is closer to z-axis and should yield higher density values
  // New maximum is therefore at z0Trk2
  BOOST_CHECK_EQUAL(*res2, z0Trk2);

  auto resWidth1 = grid.getMaxZPositionAndWidth(mainDensityMap);
  BOOST_CHECK(resWidth1.ok());
  BOOST_CHECK_EQUAL((*resWidth1).first, z0Trk2);
  BOOST_CHECK((*resWidth1).second > 0);
}

BOOST_AUTO_TEST_CASE(adaptive_gaussian_grid_density_highest_density_sum_test) {
  const int trkGridSize = 15;

  double binSize = 0.1;  // mm

  // Set up grid density with zMinMax
  AdaptiveGridTrackDensity<trkGridSize>::Config cfg(binSize);
  cfg.useHighestSumZPosition = true;

  AdaptiveGridTrackDensity<trkGridSize> grid(cfg);

  // Create some test tracks
  Covariance covMat(Covariance::Identity());

  float z0Trk1 = 0.25;
  float z0Trk2 = -10.95;
  BoundVector paramVec1;
  paramVec1 << 0.01, z0Trk1, 0, 0, 0, 0;
  BoundVector paramVec2;
  paramVec2 << 0.009, z0Trk2, 0, 0, 0, 0;

  // Create perigee surface
  std::shared_ptr<PerigeeSurface> perigeeSurface =
      Surface::makeShared<PerigeeSurface>(Vector3(0., 0., 0.));

  BoundTrackParameters params1(perigeeSurface, paramVec1, covMat);
  BoundTrackParameters params2(perigeeSurface, paramVec2, covMat);

  // Empty map
  AdaptiveGridTrackDensity<trkGridSize>::DensityMap mainDensityMap;

  // Fill grid with track densities
  auto trackDensityMap = grid.addTrack(params1, mainDensityMap);

  auto res1 = grid.getMaxZPosition(mainDensityMap);
  BOOST_CHECK(res1.ok());
  // Maximum should be at z0Trk1 position
  BOOST_CHECK_EQUAL(*res1, z0Trk1);

  // Add second track
  trackDensityMap = grid.addTrack(params2, mainDensityMap);
  auto res2 = grid.getMaxZPosition(mainDensityMap);
  BOOST_CHECK(res2.ok());
  // Trk 2 is closer to z-axis and should yield higher density values
  // New maximum is therefore at z0Trk2
  BOOST_CHECK_EQUAL(*res2, z0Trk2);

  // Add small density values around the maximum of track 1
  const float densityToAdd = 5e-4;
  mainDensityMap.at(1) += densityToAdd;
  mainDensityMap.at(3) += densityToAdd;

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
  AdaptiveGridTrackDensity<trkGridSize>::Config cfg(binSize);
  AdaptiveGridTrackDensity<trkGridSize> grid(cfg);

  // Create some test tracks
  Covariance covMat(Covariance::Identity());

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

  BoundTrackParameters params0(perigeeSurface, paramVec0, covMat);
  BoundTrackParameters params1(perigeeSurface, paramVec1, covMat);

  // Empty map
  AdaptiveGridTrackDensity<trkGridSize>::DensityMap mainDensityMap;

  // Add track 0
  auto trackDensityMap0 = grid.addTrack(params0, mainDensityMap);
  BOOST_CHECK(not mainDensityMap.empty());
  // Grid size should match trkGridSize
  BOOST_CHECK_EQUAL(mainDensityMap.size(), trkGridSize);

  // Calculate total density
  float densitySum0 = 0;
  for (auto it = mainDensityMap.begin(); it != mainDensityMap.end(); it++) {
    densitySum0 += it->second;
  }

  // Add track 0 again
  trackDensityMap0 = grid.addTrack(params0, mainDensityMap);
  BOOST_CHECK(not mainDensityMap.empty());
  // Grid size should still match trkGridSize
  BOOST_CHECK_EQUAL(mainDensityMap.size(), trkGridSize);

  // Calculate new total density
  float densitySum1 = 0;
  for (auto it = mainDensityMap.begin(); it != mainDensityMap.end(); it++) {
    densitySum1 += it->second;
  }

  BOOST_CHECK(2 * densitySum0 == densitySum1);

  // Remove track 0
  grid.subtractTrack(trackDensityMap0, mainDensityMap);

  // Calculate new total density
  float densitySum2 = 0;
  for (auto it = mainDensityMap.begin(); it != mainDensityMap.end(); it++) {
    densitySum2 += it->second;
  }

  // Density should be old one again
  BOOST_CHECK(densitySum0 == densitySum2);
  // Grid size should still match trkGridSize (removal does not touch grid size)
  BOOST_CHECK_EQUAL(mainDensityMap.size(), trkGridSize);

  // Add track 1, overlapping track 0
  auto trackDensityMap1 = grid.addTrack(params1, mainDensityMap);

  int nNonOverlappingBins = int(std::abs(z0Trk1 - z0Trk2) / binSize + 1);
  BOOST_CHECK_EQUAL(mainDensityMap.size(), trkGridSize + nNonOverlappingBins);

  float densitySum3 = 0;
  for (auto it = mainDensityMap.begin(); it != mainDensityMap.end(); it++) {
    densitySum3 += it->second;
  }

  // Remove second track 0
  grid.subtractTrack(trackDensityMap0, mainDensityMap);

  float densitySum4 = 0;
  for (auto it = mainDensityMap.begin(); it != mainDensityMap.end(); it++) {
    densitySum4 += it->second;
  }

  // Density should match differences of removed tracks
  CHECK_CLOSE_ABS(densitySum4, densitySum3 - densitySum0, 1e-5);

  // Remove track 1
  grid.subtractTrack(trackDensityMap1, mainDensityMap);

  // Size should not have changed
  BOOST_CHECK_EQUAL(mainDensityMap.size(), trkGridSize + nNonOverlappingBins);

  float densitySum5 = 0;
  for (auto it = mainDensityMap.begin(); it != mainDensityMap.end(); it++) {
    densitySum5 += it->second;
  }

  // Grid is now empty after all tracks were removed
  CHECK_CLOSE_ABS(densitySum5, 0., 1e-5);
}

}  // namespace Test
}  // namespace Acts

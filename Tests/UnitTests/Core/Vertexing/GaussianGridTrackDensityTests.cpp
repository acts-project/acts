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
#include "Acts/Definitions/Units.hpp"
#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/Surfaces/PerigeeSurface.hpp"
#include "Acts/Tests/CommonHelpers/FloatComparisons.hpp"
#include "Acts/Vertexing/GaussianGridTrackDensity.hpp"

namespace bdata = boost::unit_test::data;
using namespace Acts::UnitLiterals;

namespace Acts {
namespace Test {

using Covariance = BoundSymMatrix;

// Create a test context
GeometryContext geoContext = GeometryContext();

BOOST_AUTO_TEST_CASE(gaussian_grid_density_test) {
  // Define the size of the grids
  constexpr size_t mainGridSize = 400;
  constexpr size_t trkGridSize = 15;

  using Grid = GaussianGridTrackDensity<mainGridSize, trkGridSize>;

  double binSize = 0.1;  // mm
  double zMinMax = mainGridSize / 2 * binSize;

  // Set up grid density with zMinMax
  Grid::Config cfg(zMinMax);
  Grid grid(cfg);

  // Create some test tracks
  Covariance covMat;
  covMat << 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0,
      0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1;

  BoundVector paramVec1;
  paramVec1 << 0.01, 0.15, 0, 0, 0, 0;

  BoundVector paramVec2;
  paramVec2 << trkGridSize * binSize - 0.1, 0.15, 0, 0, 0, 0;

  BoundVector paramVec3;
  paramVec3 << trkGridSize * binSize + 0.01, 0.15, 0, 0, 0, 0;

  BoundVector paramVec3_1;
  paramVec3_1 << -(trkGridSize * binSize + 0.01), 0.15, 0, 0, 0, 0;

  BoundVector paramVec4;
  paramVec4 << 0.01, 19.95, 0, 0, 0, 0;

  BoundVector paramVec5;
  paramVec5 << 0.01, -19.95, 0, 0, 0, 0;

  BoundVector paramVec6;
  paramVec6 << 0.01, -100.0, 0, 0, 0, 0;

  BoundVector paramVec7;
  paramVec7 << 0.01, +100.0, 0, 0, 0, 0;

  // Create perigee surface
  std::shared_ptr<PerigeeSurface> perigeeSurface =
      Surface::makeShared<PerigeeSurface>(Vector3(0., 0., 0.));

  BoundTrackParameters params1(perigeeSurface, paramVec1, covMat);
  BoundTrackParameters params2(perigeeSurface, paramVec2, covMat);
  BoundTrackParameters params3(perigeeSurface, paramVec3, covMat);
  BoundTrackParameters params3_1(perigeeSurface, paramVec3_1, covMat);
  BoundTrackParameters params4(perigeeSurface, paramVec4, covMat);
  BoundTrackParameters params5(perigeeSurface, paramVec5, covMat);
  BoundTrackParameters params6(perigeeSurface, paramVec6, covMat);
  BoundTrackParameters params7(perigeeSurface, paramVec7, covMat);

  // The grid to be filled
  Grid::MainGridVector mainGrid = Grid::MainGridVector::Zero();

  // addTrack method returns the central z bin where the track density
  // grid was added and the track density grid itself for caching
  std::pair<int, Grid::TrackGridVector> binAndTrackGrid;

  // Adds tracks too far away in transverse distance
  binAndTrackGrid = grid.addTrack(params3, mainGrid);
  binAndTrackGrid = grid.addTrack(params3_1, mainGrid);
  // Adds tracks too far away in longitudinal distance
  binAndTrackGrid = grid.addTrack(params6, mainGrid);
  binAndTrackGrid = grid.addTrack(params7, mainGrid);

  // Tracks are far away from z-axis (or not in region of interest) and
  // should not have contributed to density grid
  auto zeroGrid = Grid::MainGridVector::Zero();
  BOOST_CHECK_EQUAL(mainGrid, zeroGrid);

  // Now add track 1 and 2 to grid, seperately.
  binAndTrackGrid = grid.addTrack(params1, mainGrid);
  auto gridCopy = mainGrid;

  mainGrid = Grid::MainGridVector::Zero();
  binAndTrackGrid = grid.addTrack(params2, mainGrid);

  // Track 1 is closer to z-axis and should thus yield higher
  // density values
  BOOST_CHECK(gridCopy.sum() > mainGrid.sum());

  // Track 1 and 2 summed should give higher densities than
  // only track 1 alone
  binAndTrackGrid = grid.addTrack(params1, mainGrid);
  BOOST_CHECK_EQUAL(gridCopy.sum(), mainGrid.sum());

  binAndTrackGrid = grid.addTrack(params4, mainGrid);

  // Check upper boundary
  BOOST_CHECK_EQUAL(mainGrid(mainGridSize - int((trkGridSize - 1) / 2) - 2),
                    0.);
  BOOST_CHECK(mainGrid(mainGridSize - int((trkGridSize - 1) / 2) - 1) > 0.);
  BOOST_CHECK(mainGrid(mainGridSize - 1) > 0.);

  binAndTrackGrid = grid.addTrack(params5, mainGrid);
  // Check lower boundary
  BOOST_CHECK_EQUAL(mainGrid(int((trkGridSize - 1) / 2) + 1), 0.);
  BOOST_CHECK(mainGrid(int((trkGridSize - 1) / 2)) > 0.);
  BOOST_CHECK(mainGrid(0) > 0.);

  // Check if position of maximum is correct
  auto maxRes = grid.getMaxZPosition(mainGrid);
  int maxBin = static_cast<int>((*maxRes / binSize) + mainGridSize / 2);
  BOOST_CHECK_EQUAL(maxBin, 0);

  // Check if error is thrown for empty grid
  mainGrid = Grid::MainGridVector::Zero();
  auto maxResErr = grid.getMaxZPosition(mainGrid);
  BOOST_CHECK(!maxResErr.ok());

  // Check if removal of tracks works as desired
  binAndTrackGrid = grid.addTrack(params1, mainGrid);
  binAndTrackGrid = grid.addTrack(params2, mainGrid);
  // Copy grid for future reference
  gridCopy = mainGrid;
  binAndTrackGrid = grid.addTrack(params4, mainGrid);
  // Main grid should have changed by adding track4
  BOOST_CHECK(gridCopy != mainGrid);
  // Remove track 4 again
  int zBin = binAndTrackGrid.first;
  auto trackGrid = binAndTrackGrid.second;
  grid.removeTrackGridFromMainGrid(zBin, trackGrid, mainGrid);
  // Check if it works
  BOOST_CHECK_EQUAL(gridCopy, mainGrid);
}

/// @brief Tests the functionality of the `useHighestSumZPosition` option
BOOST_AUTO_TEST_CASE(gaussian_grid_sum_max_densitytest) {
  // Define the size of the grids
  constexpr int mainGridSize = 50;
  constexpr int trkGridSize = 11;

  using Grid = Acts::GaussianGridTrackDensity<mainGridSize, trkGridSize>;

  double binSize = 0.1;  // mm
  double zMinMax = mainGridSize / 2 * binSize;

  // Set up grid density with zMinMax
  Grid::Config cfg(zMinMax);
  cfg.useHighestSumZPosition = true;
  Grid grid(cfg);

  // Create some test tracks
  Covariance covMat;
  covMat << 1e-2, 0, 0, 0, 0, 0, 0, 1e-2, 0, 0, 0, 0, 0, 0, 1e-2, 0, 0, 0, 0, 0,
      0, 1e-2, 0, 0, 0, 0, 0, 0, 1e-2, 0, 0, 0, 0, 0, 0, 1e-2;

  const double posZ1 = -1.75;
  const double posZ2 = 1.75;

  // Take two tracks, track 1 is closer in d0 and will thus have a slightly
  // higher density
  BoundVector paramVec1;
  paramVec1 << 0.01, posZ1, 0, 0, 0, 0;
  BoundVector paramVec2;
  paramVec2 << 0.015, posZ2, 0, 0, 0, 0;

  // Create perigee surface
  std::shared_ptr<PerigeeSurface> perigeeSurface =
      Surface::makeShared<PerigeeSurface>(Vector3(0., 0., 0.));

  BoundTrackParameters params1(perigeeSurface, paramVec1, covMat);
  BoundTrackParameters params2(perigeeSurface, paramVec2, covMat);

  // The grid to be filled
  Grid::MainGridVector mainGrid = Grid::MainGridVector::Zero();

  // addTrack method returns the central z bin where the track density
  // grid was added and the track density grid itself for caching
  std::pair<int, Grid::TrackGridVector> binAndTrackGrid;

  binAndTrackGrid = grid.addTrack(params1, mainGrid);
  binAndTrackGrid = grid.addTrack(params2, mainGrid);

  // Artifically add some more density around the peak of track 2
  int maxZbin = static_cast<int>((posZ2 / binSize + mainGridSize / 2.));
  mainGrid(maxZbin - 1) += 1;
  mainGrid(maxZbin + 1) += 1;

  // Even though peak density of track 1 is slighly higher, track 2
  // has a higher sum of track densities including the peak and the two
  // surrounding bins and will be the output z position.
  auto maxRes = grid.getMaxZPosition(mainGrid);
  BOOST_CHECK(maxRes.ok());
  BOOST_CHECK_EQUAL(*maxRes, posZ2);
}

/// @brief Tests the seed width
BOOST_AUTO_TEST_CASE(gaussian_grid_seed_width_test) {
  // Define the size of the grids
  constexpr int mainGridSize = 50;
  constexpr int trkGridSize = 11;

  using Grid = Acts::GaussianGridTrackDensity<mainGridSize, trkGridSize>;

  double binSize = 0.1;  // mm
  double zMinMax = mainGridSize / 2 * binSize;

  // Set up grid density with zMinMax
  Grid::Config cfg(zMinMax);
  cfg.useHighestSumZPosition = true;
  Grid grid(cfg);

  // Create some test tracks
  Covariance covMat;
  covMat << 1e-2, 0, 0, 0, 0, 0, 0, 1e-2, 0, 0, 0, 0, 0, 0, 1e-2, 0, 0, 0, 0, 0,
      0, 1e-2, 0, 0, 0, 0, 0, 0, 1e-2, 0, 0, 0, 0, 0, 0, 1e-2;

  const double posZ1 = -1.75;
  const double posZ2 = 1.75;

  // Take two tracks, track 1 is closer in d0 and will thus have a slightly
  // higher density
  BoundVector paramVec1;
  paramVec1 << 0.01, posZ1, 0, 0, 0, 0;
  BoundVector paramVec2;
  paramVec2 << 0.015, posZ2, 0, 0, 0, 0;

  // Create perigee surface
  std::shared_ptr<PerigeeSurface> perigeeSurface =
      Surface::makeShared<PerigeeSurface>(Vector3(0., 0., 0.));

  BoundTrackParameters params1(perigeeSurface, paramVec1, covMat);
  BoundTrackParameters params2(perigeeSurface, paramVec2, covMat);

  // The grid to be filled
  Grid::MainGridVector mainGrid = Grid::MainGridVector::Zero();

  // addTrack method returns the central z bin where the track density
  // grid was added and the track density grid itself for caching
  std::pair<int, Grid::TrackGridVector> binAndTrackGrid;

  binAndTrackGrid = grid.addTrack(params1, mainGrid);
  binAndTrackGrid = grid.addTrack(params2, mainGrid);

  // Artifically add some more density around the peak of track 2
  int maxZbin = static_cast<int>((posZ2 / binSize + mainGridSize / 2.));
  mainGrid(maxZbin - 1) += 1;
  mainGrid(maxZbin + 1) += 1;

  // Even though peak density of track 1 is slighly higher, track 2
  // has a higher sum of track densities including the peak and the two
  // surrounding bins and will be the output z position.

  auto maxRes = grid.getMaxZPositionAndWidth(mainGrid);
  BOOST_CHECK(maxRes.ok());
  double z = (*maxRes).first;
  double width = (*maxRes).second;

  BOOST_CHECK_EQUAL(z, posZ2);
  // Check that width was estimated
  BOOST_CHECK(width != 0.);
}

}  // namespace Test
}  // namespace Acts

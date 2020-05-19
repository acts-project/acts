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

namespace bdata = boost::unit_test::data;
using namespace Acts::UnitLiterals;

namespace Acts {
namespace Test {

using Covariance = BoundSymMatrix;

// Create a test context
GeometryContext geoContext = GeometryContext();

BOOST_AUTO_TEST_CASE(gaussian_grid_density_test) {
  // Define the size of the grids
  const int mainGridSize = 400;
  const int trkGridSize = 15;

  double binSize = 0.1;  // mm
  double zMinMax = mainGridSize / 2 * binSize;

  // Set up grid density with zMinMax
  GaussianGridTrackDensity<mainGridSize, trkGridSize>::Config cfg(zMinMax);
  GaussianGridTrackDensity<mainGridSize, trkGridSize> grid(cfg);

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
      Surface::makeShared<PerigeeSurface>(Vector3D(0., 0., 0.));

  BoundParameters params1(geoContext, covMat, paramVec1, perigeeSurface);
  BoundParameters params2(geoContext, covMat, paramVec2, perigeeSurface);
  BoundParameters params3(geoContext, covMat, paramVec3, perigeeSurface);
  BoundParameters params3_1(geoContext, covMat, paramVec3_1, perigeeSurface);
  BoundParameters params4(geoContext, covMat, paramVec4, perigeeSurface);
  BoundParameters params5(geoContext, covMat, paramVec5, perigeeSurface);
  BoundParameters params6(geoContext, covMat, paramVec6, perigeeSurface);
  BoundParameters params7(geoContext, covMat, paramVec7, perigeeSurface);

  // The grid to be filled
  ActsVectorF<mainGridSize> mainGrid(ActsVectorF<mainGridSize>::Zero());

  // addTrack method returns the central z bin where the track density
  // grid was added and the track density grid itself for caching
  std::pair<int, Acts::ActsVectorF<trkGridSize>> binAndTrackGrid;

  // Adds tracks too far away in transverse distance
  binAndTrackGrid = grid.addTrack(params3, mainGrid);
  binAndTrackGrid = grid.addTrack(params3_1, mainGrid);
  // Adds tracks too far away in longitudinal distance
  binAndTrackGrid = grid.addTrack(params6, mainGrid);
  binAndTrackGrid = grid.addTrack(params7, mainGrid);

  // Tracks are far away from z-axis (or not in region of interest) and
  // should not have contributed to density grid
  BOOST_CHECK_EQUAL(mainGrid, ActsVectorF<mainGridSize>::Zero());

  // Now add track 1 and 2 to grid, seperately.
  binAndTrackGrid = grid.addTrack(params1, mainGrid);
  auto gridCopy = mainGrid;

  mainGrid = ActsVectorF<mainGridSize>::Zero();
  binAndTrackGrid = grid.addTrack(params2, mainGrid);

  // Track 1 is closer to z-axis and should thus yield higher
  // density values
  BOOST_CHECK(gridCopy.sum() > mainGrid.sum());

  // Track 1 and 2 summed should give higher densities than
  // only track 1 alone
  binAndTrackGrid = grid.addTrack(params1, mainGrid);
  BOOST_CHECK(gridCopy.sum() < mainGrid.sum());

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
  int maxBin = (*maxRes / binSize) + mainGridSize / 2;
  BOOST_CHECK_EQUAL(maxBin, mainGridSize / 2 + 1);

  // Check if error is thrown for empty grid
  mainGrid = ActsVectorF<mainGridSize>::Zero();
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
  const int mainGridSize = 50;
  const int trkGridSize = 11;

  double binSize = 0.1;  // mm
  double zMinMax = mainGridSize / 2 * binSize;

  // Set up grid density with zMinMax
  GaussianGridTrackDensity<mainGridSize, trkGridSize>::Config cfg(zMinMax);
  cfg.useHighestSumZPosition = true;
  GaussianGridTrackDensity<mainGridSize, trkGridSize> grid(cfg);

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
      Surface::makeShared<PerigeeSurface>(Vector3D(0., 0., 0.));

  BoundParameters params1(geoContext, covMat, paramVec1, perigeeSurface);
  BoundParameters params2(geoContext, covMat, paramVec2, perigeeSurface);

  // The grid to be filled
  ActsVectorF<mainGridSize> mainGrid(ActsVectorF<mainGridSize>::Zero());

  // addTrack method returns the central z bin where the track density
  // grid was added and the track density grid itself for caching
  std::pair<int, Acts::ActsVectorF<trkGridSize>> binAndTrackGrid;

  binAndTrackGrid = grid.addTrack(params1, mainGrid);
  binAndTrackGrid = grid.addTrack(params2, mainGrid);

  // Artifically add some more density around the peak of track 2
  int maxZbin = (posZ2 / binSize + mainGridSize / 2.);
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
  const int mainGridSize = 50;
  const int trkGridSize = 11;

  double binSize = 0.1;  // mm
  double zMinMax = mainGridSize / 2 * binSize;

  // Set up grid density with zMinMax
  GaussianGridTrackDensity<mainGridSize, trkGridSize>::Config cfg(zMinMax);
  cfg.useHighestSumZPosition = true;
  GaussianGridTrackDensity<mainGridSize, trkGridSize> grid(cfg);

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
      Surface::makeShared<PerigeeSurface>(Vector3D(0., 0., 0.));

  BoundParameters params1(geoContext, covMat, paramVec1, perigeeSurface);
  BoundParameters params2(geoContext, covMat, paramVec2, perigeeSurface);

  // The grid to be filled
  ActsVectorF<mainGridSize> mainGrid(ActsVectorF<mainGridSize>::Zero());

  // addTrack method returns the central z bin where the track density
  // grid was added and the track density grid itself for caching
  std::pair<int, Acts::ActsVectorF<trkGridSize>> binAndTrackGrid;

  binAndTrackGrid = grid.addTrack(params1, mainGrid);
  binAndTrackGrid = grid.addTrack(params2, mainGrid);

  // Artifically add some more density around the peak of track 2
  int maxZbin = (posZ2 / binSize + mainGridSize / 2.);
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

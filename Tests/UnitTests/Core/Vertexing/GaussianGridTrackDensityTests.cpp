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
  Covariance covMat1;
  covMat1 << 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1,
      0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1;

  Covariance covMat2;
  covMat2 << 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1,
      0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1;

  Covariance covMat3;
  covMat3 << 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1,
      0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1;

  Covariance covMat3_1;
  covMat3_1 << 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1,
      0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1;

  Covariance covMat4;
  covMat4 << 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1,
      0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1;

  Covariance covMat5;
  covMat5 << 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1,
      0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1;

  Covariance covMat6;
  covMat6 << 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1,
      0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1;

  Covariance covMat7;
  covMat7 << 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1,
      0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1;

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

  BoundParameters params1(geoContext, std::move(covMat1), paramVec1,
                          perigeeSurface);
  BoundParameters params2(geoContext, std::move(covMat2), paramVec2,
                          perigeeSurface);
  BoundParameters params3(geoContext, std::move(covMat3), paramVec3,
                          perigeeSurface);
  BoundParameters params3_1(geoContext, std::move(covMat3_1), paramVec3_1,
                            perigeeSurface);
  BoundParameters params4(geoContext, std::move(covMat4), paramVec4,
                          perigeeSurface);
  BoundParameters params5(geoContext, std::move(covMat5), paramVec5,
                          perigeeSurface);
  BoundParameters params6(geoContext, std::move(covMat6), paramVec6,
                          perigeeSurface);
  BoundParameters params7(geoContext, std::move(covMat7), paramVec7,
                          perigeeSurface);

  // The grid to be filled
  ActsVectorF<mainGridSize> mainGrid(ActsVectorF<mainGridSize>::Zero());

  // Adds tracks too far away in transverse distance
  grid.addTrack(params3, mainGrid);
  grid.addTrack(params3_1, mainGrid);
  // Adds tracks too far away in longitudinal distance
  grid.addTrack(params6, mainGrid);
  grid.addTrack(params7, mainGrid);

  // Tracks are far away from z-axis (or not in region of interest) and
  // should not have contributed to density grid
  BOOST_CHECK_EQUAL(mainGrid, ActsVectorF<mainGridSize>::Zero());

  // Now add track 1 and 2 to grid, seperately.
  grid.addTrack(params1, mainGrid);
  auto gridCopy = mainGrid;

  mainGrid = ActsVectorF<mainGridSize>::Zero();
  grid.addTrack(params2, mainGrid);

  // Track 1 is closer to z-axis and should thus yield higher
  // density values
  BOOST_CHECK(gridCopy.sum() > mainGrid.sum());

  // Track 1 and 2 summed should give higher densities than
  // only track 1 alone
  grid.addTrack(params1, mainGrid);
  BOOST_CHECK(gridCopy.sum() < mainGrid.sum());

  grid.addTrack(params4, mainGrid);

  // Check upper boundary
  BOOST_CHECK_EQUAL(mainGrid(mainGridSize - int((trkGridSize - 1) / 2) - 2),
                    0.);
  BOOST_CHECK(mainGrid(mainGridSize - int((trkGridSize - 1) / 2) - 1) > 0.);
  BOOST_CHECK(mainGrid(mainGridSize - 1) > 0.);

  grid.addTrack(params5, mainGrid);
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
}

}  // namespace Test
}  // namespace Acts

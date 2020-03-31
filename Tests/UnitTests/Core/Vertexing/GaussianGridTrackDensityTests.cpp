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

BOOST_AUTO_TEST_CASE(testtest) {
  // Define the size of the grids
  const int mainGridSize = 400;
  const int trkGridSize = 15;

  double binSize = 0.1;  // mm
  double zMinMax = mainGridSize / 2 * binSize;

  // Set up grid density with zMinMax
  GaussianGridTrackDensity<mainGridSize, trkGridSize>::Config cfg(zMinMax);
  GaussianGridTrackDensity<mainGridSize, trkGridSize> grid(cfg);

  // Create 3 test tracks
  Covariance covMat1;
  covMat1 << 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1,
      0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1;

  Covariance covMat2;
  covMat2 << 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1,
      0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1;

  Covariance covMat3;
  covMat3 << 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1,
      0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1;

  BoundVector paramVec1;
  paramVec1 << 0.01, 0.15, 0, 0, 0, 0;

  BoundVector paramVec2;
  paramVec2 << trkGridSize * binSize - 0.1, 0.15, 0, 0, 0, 0;

  BoundVector paramVec3;
  paramVec3 << trkGridSize * binSize + 0.01, 0.15, 0, 0, 0, 0;

  // Create perigee surface
  std::shared_ptr<PerigeeSurface> perigeeSurface =
      Surface::makeShared<PerigeeSurface>(Vector3D(0., 0., 0.));

  BoundParameters params1(geoContext, std::move(covMat1), paramVec1,
                          perigeeSurface);
  BoundParameters params2(geoContext, std::move(covMat2), paramVec2,
                          perigeeSurface);
  BoundParameters params3(geoContext, std::move(covMat3), paramVec3,
                          perigeeSurface);

  // The grid to be filled
  ActsVectorF<mainGridSize> mainGrid(ActsVectorF<mainGridSize>::Zero());

  std::cout << "add track 3 ..." << std::endl;
  grid.addTrack(params3, mainGrid);

  // Track 3 is far away from z-axis and should not have contributed to
  // density grid
  BOOST_CHECK_EQUAL(mainGrid, ActsVectorF<mainGridSize>::Zero());

  std::cout << "add track 1 ..." << std::endl;
  // Now add track 1 and 2 to grid, seperately.
  grid.addTrack(params1, mainGrid);
  auto gridCopy = mainGrid;

  std::cout << "add track 2 ..." << std::endl;
  mainGrid = ActsVectorF<mainGridSize>::Zero();
  grid.addTrack(params2, mainGrid);

  // Track 1 is closer to z-axis and should thus yield higher
  // density values
  BOOST_CHECK(gridCopy.sum() > mainGrid.sum());

  // Track 1 and 2 summed should give higher densities than
  // only track 1 alone
  grid.addTrack(params1, mainGrid);
  BOOST_CHECK(gridCopy.sum() < mainGrid.sum());
}

}  // namespace Test
}  // namespace Acts

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
#include "Acts/Vertexing/AdaptiveGridTrackDensity.hpp"

namespace bdata = boost::unit_test::data;
using namespace Acts::UnitLiterals;

namespace Acts {
namespace Test {

using Covariance = BoundSymMatrix;

// Create a test context
GeometryContext geoContext = GeometryContext();

BOOST_AUTO_TEST_CASE(adaptive_gaussian_grid_density_test) {

  const int trkGridSize = 15;

  double binSize = 0.1;  // mm

  // Set up grid density with zMinMax
  AdaptiveGridTrackDensity<trkGridSize>::Config cfg(binSize);
  AdaptiveGridTrackDensity<trkGridSize> grid(cfg);

  // Create some test tracks
  Covariance covMat;
  covMat << 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0,
      0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1;

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
      Surface::makeShared<PerigeeSurface>(Vector3D(0., 0., 0.));

  BoundTrackParameters params0(perigeeSurface, paramVec0, covMat);
  BoundTrackParameters params1(perigeeSurface, paramVec1, covMat);
  BoundTrackParameters params2(perigeeSurface, paramVec2, covMat);
  BoundTrackParameters params3(perigeeSurface, paramVec3, covMat);
  BoundTrackParameters params4(perigeeSurface, paramVec4, covMat);
  BoundTrackParameters params5(perigeeSurface, paramVec5, covMat);

  // Start with empty grids
  std::vector<float> mainGridDensity;
  std::vector<int> mainGridZValues;

  // Track is too far away from z axis and was not added
  grid.addTrack(params0, mainGridDensity, mainGridZValues);
  BOOST_CHECK(mainGridDensity.empty());
  BOOST_CHECK_EQUAL(mainGridDensity.size(), mainGridZValues.size());

  // Track should have been entirely added to both grids
  grid.addTrack(params1, mainGridDensity, mainGridZValues);
  BOOST_CHECK_EQUAL(mainGridDensity.size(), trkGridSize);
  BOOST_CHECK_EQUAL(mainGridDensity.size(), mainGridZValues.size());

  // Track should have been entirely added to both grids
  grid.addTrack(params2, mainGridDensity, mainGridZValues);
  BOOST_CHECK_EQUAL(mainGridDensity.size(), 2*trkGridSize);
  BOOST_CHECK_EQUAL(mainGridDensity.size(), mainGridZValues.size());
  
  // Track 3 has overlap of 2 bins with track 1
  grid.addTrack(params3, mainGridDensity, mainGridZValues);
  BOOST_CHECK_EQUAL(mainGridDensity.size(), 3*trkGridSize - 2);
  BOOST_CHECK_EQUAL(mainGridDensity.size(), mainGridZValues.size());

  // Add first track again, should *not* introduce new z entries
  grid.addTrack(params1, mainGridDensity, mainGridZValues);
  BOOST_CHECK_EQUAL(mainGridDensity.size(), 3*trkGridSize - 2);
  BOOST_CHECK_EQUAL(mainGridDensity.size(), mainGridZValues.size());

  // Add two more tracks and check if order is correct
  grid.addTrack(params4, mainGridDensity, mainGridZValues);
  BOOST_CHECK_EQUAL(mainGridDensity.size(), mainGridZValues.size());
  grid.addTrack(params5, mainGridDensity, mainGridZValues);
  BOOST_CHECK_EQUAL(mainGridDensity.size(), mainGridZValues.size());
 
  BOOST_CHECK(std::is_sorted(std::begin(mainGridZValues), std::end(mainGridZValues)));


  // for(int i = 0; i < mainGridDensity.size(); i++){
  //   std::cout << mainGridZValues[i] << ": " << mainGridDensity[i] << std::endl;
  // }

  //std::cout << mainGridZValues << std::endl;
  // std::cout << mainGridDensity << std::endl;

}

}  // namespace Test
}  // namespace Acts

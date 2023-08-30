// This file is part of the Acts project.
//
// Copyright (C) 2020-2023 CERN for the benefit of the Acts project
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

using Covariance = BoundSquareMatrix;

BOOST_AUTO_TEST_CASE(compare_to_analytical_solution_for_single_track) {
  using Vector2 = Eigen::Matrix<float, 2, 1>;
  using Matrix2 = Eigen::Matrix<float, 2, 2>;
  // Using a large track grid so we can choose a small bin size
  const int trkGridSize = 4001;
  // Arbitrary (but small) bin size
  const float binSize = 3.1e-4;
  // Arbitrary impact parameters
  const float d0 = 0.4;
  const float z0 = -0.2;
  Vector2 impactParameters{d0, z0};

  Covariance covMat(Covariance::Identity() * 0.05);
  covMat(0, 1) = -0.02;
  covMat(1, 0) = -0.02;
  Matrix2 subCovMat = covMat.block<2, 2>(0, 0).cast<float>();
  BoundVector paramVec;
  paramVec << d0, z0, 0, 0, 0, 0;

  // Create perigee surface
  std::shared_ptr<PerigeeSurface> perigeeSurface =
      Surface::makeShared<PerigeeSurface>(Vector3(0., 0., 0.));

  BoundTrackParameters params1(perigeeSurface, paramVec, covMat);

  AdaptiveGridTrackDensity<trkGridSize>::Config cfg(binSize);
  AdaptiveGridTrackDensity<trkGridSize> grid(cfg);

  // Empty map
  AdaptiveGridTrackDensity<trkGridSize>::DensityMap mainDensityMap;

  // Add track
  auto trackDensityMap = grid.addTrack(params1, mainDensityMap);

  float relTol = 1e-5;
  float small = 1e-5;

  auto gaussian2D = [&](const Vector2& args, const Vector2& mus,
                        const Matrix2& sigmas) {
    Vector2 diffs = args - mus;
    float coef = 1 / std::sqrt(sigmas.determinant());
    float expo = -0.5 * diffs.transpose().dot(sigmas.inverse() * diffs);
    return coef * std::exp(expo);
  };

  for (auto const& it : mainDensityMap) {
    // Extract variables for better readability
    int zBin = it.first;
    float density = it.second;
    // Argument for 2D gaussian
    Vector2 dzVec{0., grid.getBinCenter(zBin)};
    // Compute correct density...
    float correctDensity = gaussian2D(dzVec, impactParameters, subCovMat);
    // ... and check if our result is equivalent
    CHECK_CLOSE_OR_SMALL(density, correctDensity, relTol, small);
  }

  // Analytical maximum of the Gaussian (can be obtained by expressing the
  // exponent as (az - b)^2 + c and noting correctMaxZ = b/a)
  float correctMaxZ =
      -0.5 * (subCovMat(0, 1) + subCovMat(1, 0)) / subCovMat(0, 0) * d0 + z0;
  // Analytical FWHM of the Gaussian (result similar to
  // https://en.wikipedia.org/wiki/Full_width_at_half_maximum#Normal_distribution
  // but the calculation needs to be slightly modified in our case)
  float correctFWHM = 2. * std::sqrt(2 * std::log(2.) *
                                     subCovMat.determinant() / subCovMat(0, 0));

  // Estimate maximum z position and seed width
  auto res = grid.getMaxZPositionAndWidth(mainDensityMap);
  BOOST_CHECK(res.ok());

  // Extract variables for better readability...
  float maxZ = res.value().first;
  float fwhm = res.value().second * 2.355f;

  // ... and check if they are correct (note: the optimization is not as exact
  // as the density values).
  float relTolOptimization = 1e-3;
  CHECK_CLOSE_REL(maxZ, correctMaxZ, relTolOptimization);
  CHECK_CLOSE_REL(fwhm, correctFWHM, relTolOptimization);
}

BOOST_AUTO_TEST_CASE(check_seed_width_estimation) {
  // Dummy track grid size (not needed for this unit test)
  const int trkGridSize = 1;
  float binSize = 2.;
  AdaptiveGridTrackDensity<trkGridSize>::Config cfg(binSize);
  AdaptiveGridTrackDensity<trkGridSize> grid(cfg);

  // Empty map
  AdaptiveGridTrackDensity<trkGridSize>::DensityMap mainDensityMap;

  // z-position of the maximum density
  float correctMaxZ = -2.;

  // Fill map with a triangular track density.
  // We use an isoscele triangle with a maximum density value of 1 and a width
  // of 20 mm. The linear approximation we use during the seed width estimation
  // should be exact in this case.
  for (int i = -6; i <= 4; i++) {
    mainDensityMap[i] =
        1.0 - 0.1 * std::abs(correctMaxZ - grid.getBinCenter(i));
  }

  // Get maximum z position and corresponding seed width
  auto res = grid.getMaxZPositionAndWidth(mainDensityMap);
  BOOST_CHECK(res.ok());

  // Check if we found the correct maximum
  float maxZ = res.value().first;
  BOOST_CHECK_EQUAL(maxZ, correctMaxZ);

  // Calculate full width at half maximum from seed width and check if it's
  // correct
  float fwhm = res.value().second * 2.355f;
  BOOST_CHECK_EQUAL(fwhm, 10.);
}

BOOST_AUTO_TEST_CASE(adaptive_gaussian_grid_density_track_adding_test) {
  const int trkGridSize = 15;

  double binSize = 0.1;  // mm

  AdaptiveGridTrackDensity<trkGridSize>::Config cfg(binSize);
  AdaptiveGridTrackDensity<trkGridSize> grid(cfg);

  // Create some test tracks in such a way that some tracks
  //  e.g. overlap and that certain tracks need to be inserted
  // between two other tracks
  Covariance covMat(Covariance::Identity());

  BoundVector paramVec0;
  paramVec0 << 100.0, -0.4, 0, 0, 0, 0;
  BoundVector paramVec1;
  paramVec1 << 0.01, -0.4, 0, 0, 0, 0;
  BoundVector paramVec2;
  paramVec2 << 0.01, 10.9, 0, 0, 0, 0;
  BoundVector paramVec3;
  paramVec3 << 0.01, 0.9, 0, 0, 0, 0;

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
  const int trkGridSize = 29;

  double binSize = 0.05;  // mm

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
  const int trkGridSize = 29;

  double binSize = 0.05;  // mm

  AdaptiveGridTrackDensity<trkGridSize>::Config cfg(binSize);
  cfg.useHighestSumZPosition = true;

  AdaptiveGridTrackDensity<trkGridSize> grid(cfg);

  // Create some test tracks
  Covariance covMat(Covariance::Identity() * 0.005);

  float z0Trk1 = 0.25;
  float z0Trk2 = -10.95;
  BoundVector paramVec1;
  paramVec1 << 0.01, z0Trk1, 0, 0, 0, 0;
  BoundVector paramVec2;
  paramVec2 << 0.0095, z0Trk2, 0, 0, 0, 0;

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
  const float densityToAdd = 0.5;
  mainDensityMap.at(4) += densityToAdd;
  mainDensityMap.at(6) += densityToAdd;

  auto res3 = grid.getMaxZPosition(mainDensityMap);
  BOOST_CHECK(res3.ok());
  // Trk 2 still has the highest peak density value, however, the small
  // added densities for track 1 around its maximum should now lead to
  // a prediction at z0Trk1 again
  BOOST_CHECK_EQUAL(*res3, z0Trk1);
}

BOOST_AUTO_TEST_CASE(adaptive_gaussian_grid_density_track_removing_test) {
  const int trkGridSize = 29;

  double binSize = 0.05;  // mm

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

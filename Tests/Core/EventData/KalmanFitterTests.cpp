// This file is part of the Acts project.
//
// Copyright (C) 2016-2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#define BOOST_TEST_MODULE KalmanFitter Tests
#include <boost/test/included/unit_test.hpp>

#include <boost/test/data/test_case.hpp>

#include "Acts/Detector/TrackingGeometry.hpp"
#include "DetectorBuild.hpp"
#include "Acts/EventData/Measurement.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/Definitions.hpp"
#include <vector>

namespace bdata = boost::unit_test::data;

namespace Acts {
namespace Test {

	using id = unsigned long int;
  ///
  /// @brief Unit test for Kalman fitter
  ///
  BOOST_AUTO_TEST_CASE(kalman_fitter_initialization)
  {
    std::shared_ptr<TrackingGeometry> detector = buildGeometry();

	std::vector<FittableMeasurement<id>> measurements;
	
	ActsSymMatrixD<2> cov2D;
	cov2D << 1. * units::_mm, 0., 0., 1. * units::_mm;
	
	Vector3D pos(-2. * units::_m, 0., 0.);
	Surface const* sur = detector->lowestTrackingVolume(pos)->associatedLayer(pos)->surfaceArray()->at(pos)[0];
	measurements.push_back(Measurement<id, eLOC_0, eLOC_1>(*sur, 0, cov2D, 0., 0.));
	
	pos = {-1. * units::_m, 0., 0.};
	sur = detector->lowestTrackingVolume(pos)->associatedLayer(pos)->surfaceArray()->at(pos)[0];
	measurements.push_back(Measurement<id, eLOC_0, eLOC_1>(*sur, 1, cov2D, 0., 0.));
	
	ActsSymMatrixD<1> cov1D;
	cov1D << 1. * units::_mm;
	
	pos = {1. * units::_m - 1. * units::_mm, 0., 0.};
	sur = detector->lowestTrackingVolume(pos)->associatedLayer(pos)->surfaceArray()->at(pos)[0];
	measurements.push_back(Measurement<id, eLOC_0>(*sur, 2, cov1D, 0.));
	
	pos = {1. * units::_m + 1. * units::_mm, 0., 0.};
	sur = detector->lowestTrackingVolume(pos)->associatedLayer(pos)->surfaceArray()->at(pos)[0];
	measurements.push_back(Measurement<id, eLOC_1>(*sur, 3, cov1D, 0.));
	
	pos = {2. * units::_m - 1. * units::_mm, 0., 0.};
	sur = detector->lowestTrackingVolume(pos)->associatedLayer(pos)->surfaceArray()->at(pos)[0];
	measurements.push_back(Measurement<id, eLOC_0>(*sur, 4, cov1D, 0.));
	
	pos = {2. * units::_m + 1. * units::_mm, 0., 0.};
	sur = detector->lowestTrackingVolume(pos)->associatedLayer(pos)->surfaceArray()->at(pos)[0];
	measurements.push_back(Measurement<id, eLOC_1>(*sur, 5, cov1D, 0.));
	
    BOOST_TEST(true);
  }
  
  BOOST_DATA_TEST_CASE(kalman_fitter_rnddata, bdata::random(
        (bdata::seed = 20,
         bdata::distribution = std::uniform_real_distribution<>(-0.5 * units::_m, 0.5 * units::_m))) ^
         bdata::random(
        (bdata::seed = 21,
         bdata::distribution = std::uniform_real_distribution<>(-0.5 * units::_m, 0.5 * units::_m))) ^
         bdata::random(
        (bdata::seed = 22,
         bdata::distribution = std::uniform_real_distribution<>(-0.5 * units::_m, 0.5 * units::_m))) ^
         bdata::random(
        (bdata::seed = 23,
         bdata::distribution = std::uniform_real_distribution<>(-0.5 * units::_m, 0.5 * units::_m))) ^
         bdata::random(
        (bdata::seed = 24,
         bdata::distribution = std::uniform_real_distribution<>(-0.5 * units::_m, 0.5 * units::_m))) ^
         bdata::random(
        (bdata::seed = 25,
         bdata::distribution = std::uniform_real_distribution<>(-0.5 * units::_m, 0.5 * units::_m))) ^
         bdata::random(
        (bdata::seed = 26,
         bdata::distribution = std::uniform_real_distribution<>(-0.5 * units::_m, 0.5 * units::_m))) ^
         bdata::random(
        (bdata::seed = 27,
         bdata::distribution = std::uniform_real_distribution<>(-0.5 * units::_m, 0.5 * units::_m))) ^
         bdata::xrange(1000), 
         x1, y1, x2, y2, x3, y3, x4, y4, index)
  {
    std::shared_ptr<TrackingGeometry> detector = buildGeometry();

	std::vector<FittableMeasurement<id>> measurements;
	
	ActsSymMatrixD<2> cov2D;
	cov2D << 1. * units::_mm, 0., 0., 1. * units::_mm;
	
	Vector3D pos(-2. * units::_m, 0., 0.);
	Surface const* sur = detector->lowestTrackingVolume(pos)->associatedLayer(pos)->surfaceArray()->at(pos)[0];
	measurements.push_back(Measurement<id, eLOC_0, eLOC_1>(*sur, index + 0, cov2D, x1, y1));
	
	pos = {-1. * units::_m, 0., 0.};
	sur = detector->lowestTrackingVolume(pos)->associatedLayer(pos)->surfaceArray()->at(pos)[0];
	measurements.push_back(Measurement<id, eLOC_0, eLOC_1>(*sur, index + 1, cov2D, x2, y2));
	
	ActsSymMatrixD<1> cov1D;
	cov1D << 1. * units::_mm;
	
	pos = {1. * units::_m - 1. * units::_mm, 0., 0.};
	sur = detector->lowestTrackingVolume(pos)->associatedLayer(pos)->surfaceArray()->at(pos)[0];
	measurements.push_back(Measurement<id, eLOC_0>(*sur, index + 2, cov1D, x3));
	
	pos = {1. * units::_m + 1. * units::_mm, 0., 0.};
	sur = detector->lowestTrackingVolume(pos)->associatedLayer(pos)->surfaceArray()->at(pos)[0];
	measurements.push_back(Measurement<id, eLOC_1>(*sur, index + 3, cov1D, y3));
	
	pos = {2. * units::_m - 1. * units::_mm, 0., 0.};
	sur = detector->lowestTrackingVolume(pos)->associatedLayer(pos)->surfaceArray()->at(pos)[0];
	measurements.push_back(Measurement<id, eLOC_0>(*sur, index + 4, cov1D, x4));
	
	pos = {2. * units::_m + 1. * units::_mm, 0., 0.};
	sur = detector->lowestTrackingVolume(pos)->associatedLayer(pos)->surfaceArray()->at(pos)[0];
	measurements.push_back(Measurement<id, eLOC_1>(*sur, index + 5, cov1D, y4));
	
    BOOST_TEST(true);
  }

}  // namespace Test
}  // namespace Acts

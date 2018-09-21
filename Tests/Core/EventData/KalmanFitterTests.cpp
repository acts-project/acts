// This file is part of the Acts project.
//
// Copyright (C) 2016-2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#define BOOST_TEST_MODULE KalmanFitter Tests
#include <boost/test/included/unit_test.hpp>

#include <vector>
#include "Acts/Detector/TrackingGeometry.hpp"
#include "Acts/EventData/Measurement.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Extrapolator/Navigator.hpp"
#include "Acts/Propagator/StraightLineStepper.hpp"
#include "Acts/Propagator/Propagator.hpp"
#include "Acts/EventData/SingleBoundTrackParameters.hpp"
#include "DetectorBuild.hpp"

namespace Acts {
namespace Test {

using id = unsigned long int;

struct SurfaceCollector
{
	std::vector<Surface*> surfaces;
	
	std::shared_ptr<TrackingGeometry> detector;
	std::vector<FittableMeasurement<id>> measurements;
	
    SurfaceCollector() = default;

    template <typename propagator_state_t>
    void
    operator()(propagator_state_t& state) const
    {
//~ std::cout << state.options.debugString << std::endl;
std::cout << state.stepping.position().x() << "\t" << state.stepping.position().y() << "\t" << state.stepping.position().z() << "\t" << state.navigation.currentSurface << std::endl;
std::cout << detector->lowestTrackingVolume(state.stepping.position())->volumeName() << "\t" << state.navigation.currentVolume->volumeName() << std::endl;
std::cout << state.navigation.navSurfaces.size() << "\t" << state.navigation.navBoundaries.size() << std::endl;
    }
};
  
  ///
  /// @brief Unit test for Kalman fitter with measurements along the x-axis
  ///
  BOOST_AUTO_TEST_CASE(kalman_fitter_initialization)
  {
	// Build detector
    std::shared_ptr<TrackingGeometry> detector = buildGeometry();

	// Construct measurements
    std::vector<FittableMeasurement<id>> measurements;

    ActsSymMatrixD<2> cov2D;
    cov2D << 1. * units::_mm, 0., 0., 1. * units::_mm;

    Vector3D       pos(-2. * units::_m, 0., 0.);
    Surface const* sur = detector->lowestTrackingVolume(pos)
                             ->associatedLayer(pos)
                             ->surfaceArray()
                             ->at(pos)[0];
std::cout << sur << "\t";
    measurements.push_back(
        Measurement<id, eLOC_0, eLOC_1>(*sur, 0, cov2D, 0., 0.));

    pos = {-1. * units::_m, 0., 0.};
    sur = detector->lowestTrackingVolume(pos)
              ->associatedLayer(pos)
              ->surfaceArray()
              ->at(pos)[0];
std::cout << sur << "\t";
    measurements.push_back(
        Measurement<id, eLOC_0, eLOC_1>(*sur, 1, cov2D, 0., 0.));

    ActsSymMatrixD<1> cov1D;
    cov1D << 1. * units::_mm;

    pos = {1. * units::_m - 1. * units::_mm, 0., 0.};
    sur = detector->lowestTrackingVolume(pos)
              ->associatedLayer(pos)
              ->surfaceArray()
              ->at(pos)[0];
std::cout << sur << "\t";
    measurements.push_back(Measurement<id, eLOC_0>(*sur, 2, cov1D, 0.));

    pos = {1. * units::_m + 1. * units::_mm, 0., 0.};
    sur = detector->lowestTrackingVolume(pos)
              ->associatedLayer(pos)
              ->surfaceArray()
              ->at(pos)[0];
std::cout << sur << "\t";
    measurements.push_back(Measurement<id, eLOC_1>(*sur, 3, cov1D, 0.));

    pos = {2. * units::_m - 1. * units::_mm, 0., 0.};
    sur = detector->lowestTrackingVolume(pos)
              ->associatedLayer(pos)
              ->surfaceArray()
              ->at(pos)[0];
std::cout << sur << "\t";
    measurements.push_back(Measurement<id, eLOC_0>(*sur, 4, cov1D, 0.));

    pos = {2. * units::_m + 1. * units::_mm, 0., 0.};
    sur = detector->lowestTrackingVolume(pos)
              ->associatedLayer(pos)
              ->surfaceArray()
              ->at(pos)[0];
std::cout << sur << std::endl;
    measurements.push_back(Measurement<id, eLOC_1>(*sur, 5, cov1D, 0.));

	// Build navigator
	Navigator navi(detector);
	navi.resolvePassive = true;
	navi.resolveMaterial = false;
	navi.resolveSensitive = false;
	
	StraightLineStepper sls;
	
	Propagator<StraightLineStepper, Navigator> prop(sls, navi);
	
    ActsSymMatrixD<5> cov;
    cov << 10 * units::_mm, 0, 0.123, 0, 0.5, 0, 10 * units::_mm, 0, 0.162, 0,
        0.123, 0, 0.1, 0, 0, 0, 0.162, 0, 0.1, 0, 0.5, 0, 0, 0,
        1. / (10 * units::_GeV);
    auto covPtr = std::make_unique<const ActsSymMatrixD<5>>(cov);
	
	Vector3D startParams(-3. * units::_m, 0., 0.), startMom(1. * units::_GeV, 0., 0);
	//~ Vector3D startParams(0., 0., 0.), startMom(1. * units::_GeV, 0., 0);
                             
	SingleCurvilinearTrackParameters<NeutralPolicy> sbtp(std::move(covPtr), startParams, startMom);

	ActionList<SurfaceCollector> aList;
	aList.get<SurfaceCollector>().detector = detector;
	aList.get<SurfaceCollector>().measurements = measurements;
		
    Propagator<StraightLineStepper, Navigator>::Options<ActionList<SurfaceCollector>> propOpts;
    propOpts.actionList = aList;
	propOpts.debug = true;
	
	pos = {2. * units::_m + 1. * units::_mm, 0., 0.};
	//~ pos = {-2. * units::_m, 0., 0.};
	sur = detector->lowestTrackingVolume(pos)
              ->associatedLayer(pos)
              ->surfaceArray()
              ->at(pos)[0];
	prop.propagate(sbtp, *sur, propOpts);

    BOOST_TEST(true);
  }
}  // namespace Test
}  // namespace Acts

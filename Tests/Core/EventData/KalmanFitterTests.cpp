// This file is part of the Acts project.
//
// Copyright (C) 2016-2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#define BOOST_TEST_MODULE KalmanFitter Tests
#include <boost/test/included/unit_test.hpp>

#include "Acts/Detector/TrackingGeometry.hpp"
#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Surfaces/RectangleBounds.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"
#include <vector>

namespace Acts {
namespace Test {

std::shared_ptr<TrackingGeometry>
buildGeometry()
{
	// Construct the rotation
	RotationMatrix3D rotation;
	double rotationAngle = M_PI * 0.5;
    Vector3D xPos(cos(rotationAngle), 0., sin(rotationAngle));
    Vector3D yPos(0., 1., 0.);
    Vector3D zPos(-sin(rotationAngle), 0., cos(rotationAngle));
    rotation.col(0) = xPos;
    rotation.col(1) = yPos;
    rotation.col(2) = zPos;
    
    // Set translation vectors
    double eps = 1. * units::_mm;
    std::vector<Vector3D> translations;
    translations.push_back({-2., 0., 0.});
    translations.push_back({-1., 0., 0.});
    translations.push_back({1. - eps, 0., 0.});
    translations.push_back({1. + eps, 0., 0.});
    translations.push_back({2. - eps, 0., 0.});
    translations.push_back({2. + eps, 0., 0.});
    
    // Boundaries of the surface
	std::shared_ptr<const RectangleBounds> rBounds(new RectangleBounds(0.5 * units::_m, 0.5 * units::_m));
	
	// Construct surfaces
	std::vector<PlaneSurface*> surfaces;
	for(size_t i = 0; i < translations.size(); i++)
	{
		Transform3D trafo(Transform3D::Identity() * rotation);
		trafo.translation() = translations[i];
		surfaces.push_back(new PlaneSurface(std::make_shared<const Transform3D>(trafo), rBounds));
	}

	
	
	return nullptr;
}

  ///
  /// @brief Unit test for Kalman fitter
  ///
  BOOST_AUTO_TEST_CASE(kalman_fitter_initialization)
  {
    // Build detector
    std::shared_ptr<TrackingGeometry> detector = buildGeometry();

    // Construct measurements
    std::map<Surface const*, std::vector<FittableMeasurement<id>>> measurements;

    ActsSymMatrixD<2> cov2D;
    cov2D << 1. * units::_mm, 0., 0., 1. * units::_mm;

    Vector3D       pos(-2. * units::_m, 0., 0.);
    Surface const* sur = detector->lowestTrackingVolume(pos)
                             ->associatedLayer(pos)
                             ->surfaceArray()
                             ->at(pos)[0];
    measurements[sur].push_back(
        Measurement<id, eLOC_0, eLOC_1>(*sur, 0, cov2D, 0., 0.));
    pos = {-1. * units::_m, 0., 0.};
    sur = detector->lowestTrackingVolume(pos)
              ->associatedLayer(pos)
              ->surfaceArray()
              ->at(pos)[0];
    measurements[sur].push_back(
        Measurement<id, eLOC_0, eLOC_1>(*sur, 1, cov2D, 0., 0.));

    ActsSymMatrixD<1> cov1D;
    cov1D << 1. * units::_mm;

    pos = {1. * units::_m - 1. * units::_mm, 0., 0.};
    sur = detector->lowestTrackingVolume(pos)
              ->associatedLayer(pos)
              ->surfaceArray()
              ->at(pos)[0];
    measurements[sur].push_back(Measurement<id, eLOC_0>(*sur, 2, cov1D, 0.));

    pos = {1. * units::_m + 1. * units::_mm, 0., 0.};
    sur = detector->lowestTrackingVolume(pos)
              ->associatedLayer(pos)
              ->surfaceArray()
              ->at(pos)[0];
    measurements[sur].push_back(Measurement<id, eLOC_1>(*sur, 3, cov1D, 0.));

    pos = {2. * units::_m - 1. * units::_mm, 0., 0.};
    sur = detector->lowestTrackingVolume(pos)
              ->associatedLayer(pos)
              ->surfaceArray()
              ->at(pos)[0];
    measurements[sur].push_back(Measurement<id, eLOC_0>(*sur, 4, cov1D, 0.));

    pos = {2. * units::_m + 1. * units::_mm, 0., 0.};
    sur = detector->lowestTrackingVolume(pos)
              ->associatedLayer(pos)
              ->surfaceArray()
              ->at(pos)[0];
    measurements[sur].push_back(Measurement<id, eLOC_1>(*sur, 5, cov1D, 0.));

    // Build navigator
    Navigator navi(detector);
    navi.resolvePassive   = false;
    navi.resolveMaterial  = false;
    navi.resolveSensitive = true;

    // Use default stepper
    StraightLineStepper sls;
    // Build navigator
    Propagator<StraightLineStepper, Navigator> prop(sls, navi);

    // Set initial parameters for the particle track
    ActsSymMatrixD<5> cov;
    cov << 10 * units::_mm, 0, 0.123, 0, 0.5, 0, 10 * units::_mm, 0, 0.162, 0,
        0.123, 0, 0.1, 0, 0, 0, 0.162, 0, 0.1, 0, 0.5, 0, 0, 0,
        1. / (10 * units::_GeV);
    auto     covPtr = std::make_unique<const ActsSymMatrixD<5>>(cov);
    Vector3D startParams(-3. * units::_m, 0., 0.),
        //~ Vector3D startParams(-1.0006 * units::_m, 0., 499.8 * units::_mm),
        startMom(1. * units::_GeV, 0., 0);

    SingleCurvilinearTrackParameters<NeutralPolicy> sbtp(
        std::move(covPtr), startParams, startMom);

    // Create action list for surface collection
    ActionList<SurfaceCollection, SurfaceCollector<SelectSurfaceWithHit>> aList;
    aList.get<SurfaceCollection>().measurements = measurements;
    aList.get<SurfaceCollector<SelectSurfaceWithHit>>().selector.measurements
        = measurements;
    // Set options for propagator
    Propagator<StraightLineStepper, Navigator>::
        Options<ActionList<SurfaceCollection,
                           SurfaceCollector<SelectSurfaceWithHit>>>
            propOpts;
    propOpts.actionList = aList;

    // Launch and collect
    const auto&                        result = prop.propagate(sbtp, propOpts);
    const std::vector<Surface const*>& surResult
        = result.get<typename SurfaceCollection::result_type>();
    const SurfaceCollector<SelectSurfaceWithHit>::this_result& surResult2
        = result.get<
            typename SurfaceCollector<SelectSurfaceWithHit>::result_type>();

    // Test if results match the number of measurements
    BOOST_TEST(surResult.size() == 6);
    BOOST_TEST(surResult2.collected.size() == 6);

    // Re-configure propagation with B-field
    ConstantBField               bField(Vector3D(0., 0.5 * units::_T, 0.));
    EigenStepper<ConstantBField> es(bField);
    Propagator<EigenStepper<ConstantBField>, Navigator> propB(es, navi);
    covPtr = std::make_unique<const ActsSymMatrixD<5>>(cov);
    SingleCurvilinearTrackParameters<ChargedPolicy> sbtpB(
        std::move(covPtr), startParams, startMom, 1.);
    AbortList<EndOfWorld> abortList;
    Propagator<EigenStepper<ConstantBField>, Navigator>::
        Options<ActionList<SurfaceCollection,
                           SurfaceCollector<SelectSurfaceWithHit>>,
                AbortList<EndOfWorld>>
            propOptsB;
    propOptsB.actionList     = aList;
    propOptsB.stopConditions = abortList;
    propOptsB.maxSteps       = 1e6;

    const auto& resultB = propB.propagate(sbtpB, propOptsB);
    const std::vector<Surface const*>& surResultB
        = resultB.get<typename SurfaceCollection::result_type>();
    const SurfaceCollector<SelectSurfaceWithHit>::this_result& surResultB2
        = resultB.get<
            typename SurfaceCollector<SelectSurfaceWithHit>::result_type>();

    for (size_t i = 0; i < surResultB2.collected.size(); i++)
      std::cout << surResultB2.collected[i].position.x() << " "
                << surResultB2.collected[i].position.y() << " "
                << surResultB2.collected[i].position.z() << std::endl
                << std::endl;
    BOOST_TEST(surResultB.size() == 2);
    BOOST_TEST(surResultB2.collected.size() == 2);
  }

}  // namespace Test
}  // namespace Acts

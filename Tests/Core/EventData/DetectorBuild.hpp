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
#include "Acts/Surfaces/SurfaceArray.hpp"
#include "Acts/Layers/PlaneLayer.hpp"
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
    double eps = 2. * units::_mm;
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
	std::array<PlaneSurface*, 6> surfaces;
	unsigned int i;
	for(i = 0; i < translations.size(); i++)
	{
		Transform3D trafo(Transform3D::Identity() * rotation);
		trafo.translation() = translations[i];
		surfaces[i] = new PlaneSurface(std::make_shared<const Transform3D>(trafo), rBounds);
	}
    
	std::array<LayerPtr, 6> layers;
	for(i = 0; i < 6; i++)
	{
		Transform3D trafo(Transform3D::Identity() * rotation);
		trafo.translation() = translations[i];
		
		std::unique_ptr<SurfaceArray> surArray(new SurfaceArray(surfaces[i]));
		
		layers[i] = PlaneLayer::create(std::make_shared<const Transform3D>(trafo), rBounds, std::move(surArray), 0.5 * units::_mm);
	}
	
	
	
	return nullptr;
}

}  // namespace Test
}  // namespace Acts

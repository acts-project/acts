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
	BOOST_TEST(true);
  }

}  // namespace Test
}  // namespace Acts

// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "ActsFatras/Digitization/DigitizationData.hpp"
#include "ActsFatras/Digitization/DigitizationError.hpp"
#include "ActsFatras/Digitization/ReadoutProjector.hpp"
#include "ActsFatras/EventData/Hit.hpp"
#include <Acts/Geometry/GeometryContext.hpp>
#include <Acts/Surfaces/PlaneSurface.hpp>
#include <Acts/Surfaces/RectangleBounds.hpp>
#include <Acts/Tests/CommonHelpers/FloatComparisons.hpp>
#include <Acts/Utilities/Definitions.hpp>
#include <Acts/Utilities/Helpers.hpp>
#include <Acts/Utilities/ParameterDefinitions.hpp>

namespace ActsFatras {

BOOST_AUTO_TEST_SUITE(Digitization)

BOOST_AUTO_TEST_CASE(ReadoutProjectorTest) {
  Acts::GeometryContext geoCtx;

  auto rectangleBounds = std::make_shared<Acts::RectangleBounds>(2., 3.5);
  auto planeSurface = Acts::Surface::makeShared<Acts::PlaneSurface>(
      Acts::Transform3D::Identity(), rectangleBounds);
  double thickness = 0.15;
  Acts::Vector3D normal(0., 0., 1.);

  Hit fHit(0, 1, {0., 0., 0., 0.}, {0.5, 0.5, 3.0, 0.}, {0.5, 0.5, 3.0, 0.});

  DigitizationInput dInput(fHit, geoCtx, planeSurface.get());

  ReadoutProjector projector;

  // Emulate a 3D-projection
  auto readout3d = projector.project(dInput, {0., 0., 0.}, thickness);
  BOOST_CHECK(readout3d.ok());
  auto projected3d = readout3d.value();
  CHECK_CLOSE_ABS(projected3d[0].second, 0.075, Acts::s_epsilon);
  CHECK_CLOSE_ABS(projected3d[1].second, -0.075, Acts::s_epsilon);

  // Emulate a readout at the entry (bottom) surface
  auto driftE = Acts::Vector3D(0., 0.2, -0.8).normalized();
  auto readoutE = projector.project(dInput, driftE, thickness);
  BOOST_CHECK(readoutE.ok());
  auto projectedE = readoutE.value();
  CHECK_CLOSE_ABS(projectedE[0].second, 0., Acts::s_epsilon);
  CHECK_CLOSE_ABS(projectedE[1].second, -1 * thickness / normal.dot(driftE),
                  Acts::s_epsilon);

  // Emulate a readout at the exit (top) surface
  auto driftT = Acts::Vector3D(0.3, 0.2, 0.8).normalized();
  auto readoutT = projector.project(dInput, driftT, thickness);
  BOOST_CHECK(readoutT.ok());
  auto projectedT = readoutT.value();
  CHECK_CLOSE_ABS(projectedT[0].second, 1 * thickness / normal.dot(driftT),
                  Acts::s_epsilon);
  CHECK_CLOSE_ABS(projectedT[1].second, 0., Acts::s_epsilon);
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace ActsFatras
// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "ActsFatras/Digitization/Channelizer.hpp"
#include "ActsFatras/Digitization/DigitizationData.hpp"
#include "ActsFatras/Digitization/DigitizationError.hpp"
#include "ActsFatras/EventData/Hit.hpp"
#include <Acts/Geometry/GeometryContext.hpp>
#include <Acts/Surfaces/PlaneSurface.hpp>
#include <Acts/Surfaces/RectangleBounds.hpp>
#include <Acts/Tests/CommonHelpers/FloatComparisons.hpp>
#include <Acts/Utilities/BinUtility.hpp>
#include <Acts/Utilities/Definitions.hpp>
#include <Acts/Utilities/Helpers.hpp>
#include <Acts/Utilities/ParameterDefinitions.hpp>

namespace ActsFatras {

BOOST_AUTO_TEST_SUITE(Digitization)

BOOST_AUTO_TEST_CASE(ChannelizerCartesian) {
  Acts::GeometryContext geoCtx;

  auto rectangleBounds = std::make_shared<Acts::RectangleBounds>(1., 1.);
  auto planeSurface = Acts::Surface::makeShared<Acts::PlaneSurface>(
      Acts::Transform3D::Identity(), rectangleBounds);

  // The segementation
  Acts::BinUtility pixelated(20, -1., 1., Acts::open, Acts::binX);
  pixelated += Acts::BinUtility(20, -1., 1., Acts::open, Acts::binY);

  Channelizer cl;

  // Test: Normal hit into the surface
  Hit nHit(0, 1, {0.37, 0.76, 0., 0.}, {0., 0., 1.0, 0.}, {0., 0., 1.0, 0.});
  DigitizationInput nInput(nHit, geoCtx, planeSurface.get(), pixelated);
  auto nPosition = nHit.position4().segment<2>(Acts::ePos0);
  auto nSegments = cl.segments(nInput, nPosition, nPosition);
  BOOST_CHECK(nSegments.size() == 1);
  BOOST_CHECK(nSegments[0].bin[0] == 13);
  BOOST_CHECK(nSegments[0].bin[1] == 17);

  // Test: Inclined hit into the surface - negative x direction
  Acts::Vector2D ixPositionS(0.37, 0.76);
  Acts::Vector2D ixPositionE(0.02, 0.73);
  auto ixSegments = cl.segments(nInput, ixPositionS, ixPositionE);
  BOOST_CHECK(ixSegments.size() == 4);

  // Test: Inclined hit into the surface - positive y direction
  Acts::Vector2D iyPositionS(0.37, 0.76);
  Acts::Vector2D iyPositionE(0.39, 0.91);
  auto iySegments = cl.segments(nInput, iyPositionS, iyPositionE);
  BOOST_CHECK(iySegments.size() == 3);

  // Test: Inclined hit into the surface - x/y direction
  Acts::Vector2D ixyPositionS(-0.27, 0.76);
  Acts::Vector2D ixyPositionE(-0.02, -0.73);
  auto ixySegments = cl.segments(nInput, ixyPositionS, ixyPositionE);
  BOOST_CHECK(ixySegments.size() == 18);
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace ActsFatras
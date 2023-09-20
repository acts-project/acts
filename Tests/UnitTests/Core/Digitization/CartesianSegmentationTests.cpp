// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/test/data/test_case.hpp>
#include <boost/test/unit_test.hpp>

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/Units.hpp"
#include "Acts/Digitization/CartesianSegmentation.hpp"
#include "Acts/Digitization/Segmentation.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/Surfaces/RectangleBounds.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Tests/CommonHelpers/FloatComparisons.hpp"

#include <cmath>
#include <cstdlib>
#include <memory>

namespace bdata = boost::unit_test::data;
namespace tt = boost::test_tools;
using namespace Acts::UnitLiterals;

namespace Acts {
namespace Test {

size_t nbinsx = 100;
size_t nbinsy = 200;
double hThickness = 75_um;
double lAngle = 0.1;

// Module bounds
auto moduleBounds = std::make_shared<const RectangleBounds>(5_mm, 10_mm);
CartesianSegmentation cSegmentation(moduleBounds, nbinsx, nbinsy);

// Create a test context
GeometryContext tgContext = GeometryContext();

/// @brief Unit test for the Cartesian segmentation
///
BOOST_AUTO_TEST_CASE(cartesian_segmentation) {
  // The surface vectors: positive readout, zero lorentz angle
  // -> PZL
  SurfacePtrVector boundariesPZL;
  SurfacePtrVector segSurfacesXPZL;
  SurfacePtrVector segSurfacesYPZL;

  cSegmentation.createSegmentationSurfaces(boundariesPZL, segSurfacesXPZL,
                                           segSurfacesYPZL, hThickness, 1, 0.);

  BOOST_CHECK_EQUAL(boundariesPZL.size(), 6u);

  // There's one less because of the boundary and lorentz plane
  BOOST_CHECK_EQUAL(segSurfacesXPZL.size(), size_t(nbinsx - 1));
  BOOST_CHECK_EQUAL(segSurfacesYPZL.size(), size_t(nbinsy - 1));

  // Check the boundary surfaces are thickness away
  auto centerReadoutPZL = boundariesPZL[0]->center(tgContext);
  auto centerCounterPZL = boundariesPZL[1]->center(tgContext);
  auto centerDiffPZL = centerReadoutPZL - centerCounterPZL;
  double thicknessPZL = centerDiffPZL.norm();

  CHECK_CLOSE_REL(thicknessPZL, 2 * hThickness, 10e-6);

  // The surface vectors: negative readout, zero lorentz angle
  // -> NZL
  SurfacePtrVector boundariesNZL;
  SurfacePtrVector segSurfacesXNZL;
  SurfacePtrVector segSurfacesYNZL;

  cSegmentation.createSegmentationSurfaces(boundariesNZL, segSurfacesXNZL,
                                           segSurfacesYNZL, hThickness, -1, 0.);

  BOOST_CHECK_EQUAL(boundariesNZL.size(), 6u);

  // There's one less because of the boundary and lorentz plane
  BOOST_CHECK_EQUAL(segSurfacesXNZL.size(), size_t(nbinsx - 1));
  BOOST_CHECK_EQUAL(segSurfacesYNZL.size(), size_t(nbinsy - 1));

  // Check the boundary surfaces are thickness away
  auto centerReadoutNZL = boundariesNZL[0]->center(tgContext);
  auto centerCounterNZL = boundariesNZL[1]->center(tgContext);
  auto centerDiffNZL = centerReadoutNZL - centerCounterNZL;
  double thicknessNZL = centerDiffNZL.norm();

  CHECK_CLOSE_REL(thicknessNZL, 2 * hThickness, 10e-6);

  // Check that the readout / counter surfaces are reversed
  CHECK_CLOSE_OR_SMALL(centerReadoutPZL, centerCounterNZL, 10e-6, 10e-9);
  CHECK_CLOSE_OR_SMALL(centerReadoutNZL, centerCounterPZL, 10e-6, 10e-9);

  // The surface vectors: positive readout, lorentz angle
  // -> PL
  SurfacePtrVector boundariesPL;
  SurfacePtrVector segSurfacesXPL;
  SurfacePtrVector segSurfacesYPL;

  cSegmentation.createSegmentationSurfaces(
      boundariesPL, segSurfacesXPL, segSurfacesYPL, hThickness, 1, lAngle);

  BOOST_CHECK_EQUAL(boundariesPL.size(), 6u);

  // There's one less because of the boundary and lorentz plane
  BOOST_CHECK_EQUAL(segSurfacesXPL.size(), size_t(nbinsx - 1));
  BOOST_CHECK_EQUAL(segSurfacesYPL.size(), size_t(nbinsy - 1));

  // Check the boundary surfaces are thickness away
  auto centerReadoutPL = boundariesPL[0]->center(tgContext);
  auto centerCounterPL = boundariesPL[1]->center(tgContext);
  double thicknessPL = abs((centerReadoutPL - centerCounterPL).z());

  CHECK_CLOSE_REL(thicknessPL, 2 * hThickness, 10e-6);

  // check the lorentz angle - let's take the second one
  const auto* pSurface =
      dynamic_cast<const Acts::PlaneSurface*>(segSurfacesXPL[2].get());
  BOOST_REQUIRE(pSurface != nullptr);
  auto nLorentzPlane = pSurface->normal(tgContext);

  Vector3 nNominal(1., 0., 0.);
  double tAngle = acos(nLorentzPlane.dot(nNominal));

  CHECK_CLOSE_REL(tAngle, lAngle, 0.001);
}

}  // namespace Test
}  // namespace Acts

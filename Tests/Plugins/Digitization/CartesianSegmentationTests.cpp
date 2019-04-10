// This file is part of the Acts project.
//
// Copyright (C) 2019 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

// clang-format off
#define BOOST_TEST_MODULE CartesianSegmentation Tests
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#include <boost/test/data/test_case.hpp>
// clang-format on

#include "Acts/Utilities/Units.hpp"
#include "Acts/Utilities/GeometryContext.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Surfaces/PlanarBounds.hpp"
#include "Acts/Surfaces/RectangleBounds.hpp"
#include "Acts/Plugins/Digitization/CartesianSegmentation.hpp"
#include "Acts/Tests/CommonHelpers/FloatComparisons.hpp"

namespace bdata = boost::unit_test::data;
namespace tt    = boost::test_tools;

namespace Acts {

namespace Test {

  size_t nbinsx     = 100;
  size_t nbinsy     = 200;
  double hThickness = 75 * units::_um;
  double lAngle     = 0.1;

  // Module bounds
  auto moduleBounds = std::make_shared<const RectangleBounds>(5 * units::_mm,
                                                              10 * units::_mm);
  CartesianSegmentation cSegmentation(moduleBounds, nbinsx, nbinsy);

  // Create a test context
  GeometryContext tgContext = GeometryContext();

  /// @brief Unit test for the Cartesian segmentation
  ///
  BOOST_AUTO_TEST_CASE(cartesian_segmentation)
  {

    // The surface vectors: positive readout, zero lorentz angle
    // -> PZL
    SurfacePtrVector boundariesPZL;
    SurfacePtrVector segSurfacesXPZL;
    SurfacePtrVector segSurfacesYPZL;

    cSegmentation.createSegmentationSurfaces(
        boundariesPZL, segSurfacesXPZL, segSurfacesYPZL, hThickness, 1, 0.);

    BOOST_CHECK_EQUAL(boundariesPZL.size(), 6);

    // There's one less because of the boundary and lorentz plane
    BOOST_CHECK_EQUAL(segSurfacesXPZL.size(), size_t(nbinsx - 1));
    BOOST_CHECK_EQUAL(segSurfacesYPZL.size(), size_t(nbinsy - 1));

    // Check the boundary surfaces are thickness away
    auto   centerReadoutPZL = boundariesPZL[0]->center(tgContext);
    auto   centerCounterPZL = boundariesPZL[1]->center(tgContext);
    auto   centerDiffPZL    = centerReadoutPZL - centerCounterPZL;
    double thicknessPZL     = centerDiffPZL.norm();

    CHECK_CLOSE_REL(thicknessPZL, 2 * hThickness, 10e-6);

    // The surface vectors: negative readout, zero lorentz angle
    // -> NZL
    SurfacePtrVector boundariesNZL;
    SurfacePtrVector segSurfacesXNZL;
    SurfacePtrVector segSurfacesYNZL;

    cSegmentation.createSegmentationSurfaces(
        boundariesNZL, segSurfacesXNZL, segSurfacesYNZL, hThickness, -1, 0.);

    BOOST_CHECK_EQUAL(boundariesNZL.size(), 6);

    // There's one less because of the boundary and lorentz plane
    BOOST_CHECK_EQUAL(segSurfacesXNZL.size(), size_t(nbinsx - 1));
    BOOST_CHECK_EQUAL(segSurfacesYNZL.size(), size_t(nbinsy - 1));

    // Check the boundary surfaces are thickness away
    auto   centerReadoutNZL = boundariesNZL[0]->center(tgContext);
    auto   centerCounterNZL = boundariesNZL[1]->center(tgContext);
    auto   centerDiffNZL    = centerReadoutNZL - centerCounterNZL;
    double thicknessNZL     = centerDiffNZL.norm();

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

    BOOST_CHECK_EQUAL(boundariesPL.size(), 6);

    // There's one less because of the boundary and lorentz plane
    BOOST_CHECK_EQUAL(segSurfacesXPL.size(), size_t(nbinsx - 1));
    BOOST_CHECK_EQUAL(segSurfacesYPL.size(), size_t(nbinsy - 1));

    // Check the boundary surfaces are thickness away
    auto   centerReadoutPL = boundariesPL[0]->center(tgContext);
    auto   centerCoutnerPL = boundariesPL[1]->center(tgContext);
    double thicknessPL     = abs((centerReadoutPL - centerCoutnerPL).z());

    CHECK_CLOSE_REL(thicknessPL, 2 * hThickness, 10e-6);

    // check the lorentz angle - let's take the second one
    auto nLorentzPlane = segSurfacesXPL[2]->normal(tgContext);

    Vector3D nNominal(1., 0., 0.);
    double   tAngle = acos(nLorentzPlane.dot(nNominal));

    CHECK_CLOSE_REL(tAngle, lAngle, 0.001);
  }

}  // namespace
}  // namespace
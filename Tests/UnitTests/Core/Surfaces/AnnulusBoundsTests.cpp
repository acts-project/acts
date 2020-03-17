// This file is part of the Acts project.
//
// Copyright (C) 2017-2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

// clang-format off
#include <boost/test/unit_test.hpp>
#include <boost/test/data/test_case.hpp>
#include <boost/test/output_test_stream.hpp>
// clang-format on

#include <limits>
#include <fstream>

#include "Acts/Surfaces/AnnulusBounds.hpp"
#include "Acts/Tests/CommonHelpers/FloatComparisons.hpp"
#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Utilities/Helpers.hpp"

namespace Acts {

namespace Test {
BOOST_AUTO_TEST_SUITE(Surfaces)

double minRadius = 7.2;
double maxRadius = 12.0;
double minPhi = 0.74195;
double maxPhi = 1.33970;

Vector2D offset(-2., 2.);

/// Unit tests for AnnulusBounds constrcuctors
BOOST_AUTO_TEST_CASE(AnnulusBoundsConstruction) {
  //
  /// Test construction with radii and default sector
  auto original = AnnulusBounds(minRadius, maxRadius, minPhi, maxPhi, offset);
  BOOST_CHECK_EQUAL(original.type(), SurfaceBounds::Annulus);

  AnnulusBounds copied(original);
  BOOST_CHECK_EQUAL(copied.type(), SurfaceBounds::Annulus);
}

/// Unit tests for AnnulusBounds properties
BOOST_AUTO_TEST_CASE(AnnulusBoundsProperties) {
  /// Test construction with radii and default sector
  AnnulusBounds aBounds(minRadius, maxRadius, minPhi, maxPhi, offset);

  /// Test clone
  auto pClonedAnnulusBounds = aBounds.clone();
  BOOST_CHECK_NE(pClonedAnnulusBounds, nullptr);
  delete pClonedAnnulusBounds;
  //
  /// Test type() (redundant; already used in constructor confirmation)
  BOOST_CHECK_EQUAL(aBounds.type(), SurfaceBounds::Annulus);

  /// Test positions inside/outside
  // - start from cartesian (from test drawing)
  Vector2D inSurfaceXY(7., 7.);
  Vector2D outsideXY1(5., 5.);
  Vector2D outsideXY2(10., 3.);
  Vector2D outsideXY3(10., 10.);
  Vector2D outsideXY4(4., 10.);
  std::vector<Vector2D> testPoints = {inSurfaceXY, outsideXY1, outsideXY2,
                                      outsideXY3, outsideXY4};

  auto toStripFrame = [&](const Vector2D& xy) -> Vector2D {
    auto shifted = xy + offset;
    double r = VectorHelpers::perp(shifted);
    double phi = VectorHelpers::phi(shifted);
    return Vector2D(r, phi);
  };

  BOOST_CHECK(aBounds.inside(toStripFrame(inSurfaceXY), BoundaryCheck(true)));
  BOOST_CHECK(!aBounds.inside(toStripFrame(outsideXY1), BoundaryCheck(true)));
  BOOST_CHECK(!aBounds.inside(toStripFrame(outsideXY2), BoundaryCheck(true)));
  BOOST_CHECK(!aBounds.inside(toStripFrame(outsideXY3), BoundaryCheck(true)));
  BOOST_CHECK(!aBounds.inside(toStripFrame(outsideXY4), BoundaryCheck(true)));

  /// Check radial inside
  BOOST_CHECK(!aBounds.insideRadialBounds(0.5));
  BOOST_CHECK(aBounds.insideRadialBounds(9.));
  BOOST_CHECK(!aBounds.insideRadialBounds(18.));

  /// Test rMin
  BOOST_CHECK_EQUAL(aBounds.rMin(), minRadius);
  //
  /// Test rMax
  BOOST_CHECK_EQUAL(aBounds.rMax(), maxRadius);
  /// Test phiMin
  BOOST_CHECK_EQUAL(aBounds.rMin(), minRadius);
  //
  /// Test phiMax
  BOOST_CHECK_EQUAL(aBounds.rMax(), maxRadius);
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace Test

}  // namespace Acts

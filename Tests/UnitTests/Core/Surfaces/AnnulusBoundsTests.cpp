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
#include <boost/test/tools/output_test_stream.hpp>
// clang-format on

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Surfaces/AnnulusBounds.hpp"
#include "Acts/Tests/CommonHelpers/FloatComparisons.hpp"
#include "Acts/Utilities/Helpers.hpp"

#include <fstream>
#include <limits>

namespace Acts {

namespace Test {
BOOST_AUTO_TEST_SUITE(Surfaces)

double minRadius = 7.2;
double maxRadius = 12.0;
double minPhi = 0.74195;
double maxPhi = 1.33970;

Vector2 offset(-2., 2.);

// Unit tests for AnnulusBounds constrcuctors
BOOST_AUTO_TEST_CASE(AnnulusBoundsConstruction) {
  // Test construction with radii and default sector
  auto original = AnnulusBounds(minRadius, maxRadius, minPhi, maxPhi, offset);
  AnnulusBounds copied(original);
  BOOST_CHECK_EQUAL(original, copied);
}

// Unit tests for AnnulusBounds recreation
BOOST_AUTO_TEST_CASE(AnnulusBoundsRecreation) {
  // Test construction with radii and default sector
  auto original = AnnulusBounds(minRadius, maxRadius, minPhi, maxPhi, offset);
  auto valvector = original.values();
  std::array<double, AnnulusBounds::eSize> values{};
  std::copy_n(valvector.begin(), AnnulusBounds::eSize, values.begin());
  AnnulusBounds recreated(values);
  BOOST_CHECK_EQUAL(original, recreated);
}

// Unit tests for AnnulusBounds exception throwing
BOOST_AUTO_TEST_CASE(AnnulusBoundsExcpetion) {
  // Exception for negative inner radius
  BOOST_CHECK_THROW(AnnulusBounds(-1., maxRadius, minPhi, maxPhi, offset),
                    std::logic_error);
  // Exception for negative outer radius
  BOOST_CHECK_THROW(AnnulusBounds(minRadius, -1., minPhi, maxPhi, offset),
                    std::logic_error);
  // Exception for swapped radii
  BOOST_CHECK_THROW(AnnulusBounds(maxRadius, minRadius, minPhi, maxPhi, offset),
                    std::logic_error);
  // Exception for out of range  min phi
  BOOST_CHECK_THROW(AnnulusBounds(minRadius, maxRadius, -4., maxPhi, offset),
                    std::logic_error);
  // Exception for out of range  max phi
  BOOST_CHECK_THROW(AnnulusBounds(minRadius, maxRadius, minPhi, 4., offset),
                    std::logic_error);
  // Exception for out of range  max phi
  BOOST_CHECK_THROW(AnnulusBounds(minRadius, maxRadius, maxPhi, minPhi, offset),
                    std::logic_error);
}

/// Unit tests for AnnulusBounds properties
BOOST_AUTO_TEST_CASE(AnnulusBoundsProperties) {
  /// Test construction with radii and default sector
  AnnulusBounds aBounds(minRadius, maxRadius, minPhi, maxPhi, offset);

  //
  /// Test type() (redundant; already used in constructor confirmation)
  BOOST_CHECK_EQUAL(aBounds.type(), SurfaceBounds::eAnnulus);

  /// Test positions inside/outside
  // - start from cartesian (from test drawing)
  Vector2 inSurfaceXY(7., 7.);
  Vector2 outsideXY1(5., 5.);
  Vector2 outsideXY2(10., 3.);
  Vector2 outsideXY3(10., 10.);
  Vector2 outsideXY4(4., 10.);
  std::vector<Vector2> testPoints = {inSurfaceXY, outsideXY1, outsideXY2,
                                     outsideXY3, outsideXY4};

  auto toStripFrame = [&](const Vector2& xy) -> Vector2 {
    auto shifted = xy + offset;
    double r = VectorHelpers::perp(shifted);
    double phi = VectorHelpers::phi(shifted);
    return Vector2(r, phi);
  };

  BOOST_CHECK(aBounds.inside(toStripFrame(inSurfaceXY), BoundaryCheck(true)));
  BOOST_CHECK(!aBounds.inside(toStripFrame(outsideXY1), BoundaryCheck(true)));
  BOOST_CHECK(!aBounds.inside(toStripFrame(outsideXY2), BoundaryCheck(true)));
  BOOST_CHECK(!aBounds.inside(toStripFrame(outsideXY3), BoundaryCheck(true)));
  BOOST_CHECK(!aBounds.inside(toStripFrame(outsideXY4), BoundaryCheck(true)));

  // Check radial inside
  BOOST_CHECK(!aBounds.insideRadialBounds(0.5));
  BOOST_CHECK(aBounds.insideRadialBounds(9.));
  BOOST_CHECK(!aBounds.insideRadialBounds(18.));

  // Test rMin
  BOOST_CHECK_EQUAL(aBounds.get(AnnulusBounds::eMinR), minRadius);
  // Test rMax
  BOOST_CHECK_EQUAL(aBounds.get(AnnulusBounds::eMaxR), maxRadius);
  // Test phiMin
  BOOST_CHECK_EQUAL(aBounds.get(AnnulusBounds::eMinPhiRel), minPhi);
  // Test phiMax
  BOOST_CHECK_EQUAL(aBounds.get(AnnulusBounds::eMaxPhiRel), maxPhi);
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace Test

}  // namespace Acts

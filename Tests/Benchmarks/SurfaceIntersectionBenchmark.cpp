// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/test/data/test_case.hpp>
#include <boost/test/tools/output_test_stream.hpp>
#include <boost/test/unit_test.hpp>

#include "Acts/Definitions/Units.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Surfaces/CylinderBounds.hpp"
#include "Acts/Surfaces/CylinderSurface.hpp"
#include "Acts/Surfaces/DiscSurface.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/Surfaces/RadialBounds.hpp"
#include "Acts/Surfaces/RectangleBounds.hpp"
#include "Acts/Surfaces/StrawSurface.hpp"
#include "Acts/Tests/CommonHelpers/BenchmarkTools.hpp"

#include <cmath>

namespace bdata = boost::unit_test::data;
namespace tt = boost::test_tools;
using namespace Acts::UnitLiterals;

namespace Acts {
namespace Test {

// Some randomness & number crunching
unsigned int ntests = 10;
unsigned int nrepts = 2000;
const bool boundaryCheck = false;
const bool testPlane = true;
const bool testDisc = true;
const bool testCylinder = true;
const bool testStraw = true;

// Create a test context
GeometryContext tgContext = GeometryContext();

// Create a test plane in 10 m distance
// Some random transform
Transform3 at = Transform3::Identity() * Translation3(0_m, 0_m, 10_m) *
                AngleAxis3(0.15, Vector3(1.2, 1.2, 0.12).normalized());

// Define the Plane surface
auto rb = std::make_shared<RectangleBounds>(1_m, 1_m);
auto aPlane = Surface::makeShared<PlaneSurface>(at, std::move(rb));

// Define the Disc surface
auto db = std::make_shared<RadialBounds>(0.2_m, 1.2_m);
auto aDisc = Surface::makeShared<DiscSurface>(at, std::move(db));

// Define a Cylinder surface
auto cb = std::make_shared<CylinderBounds>(10_m, 100_m);
auto aCylinder = Surface::makeShared<CylinderSurface>(at, std::move(cb));

// Define a Straw surface
auto aStraw = Surface::makeShared<StrawSurface>(at, 50_cm, 2_m);

// The origin of our attempts for plane, disc and cylinder
Vector3 origin(0., 0., 0.);

// The origin for straw/line attempts
Vector3 originStraw(0.3_m, -0.2_m, 11_m);

template <typename surface_t>
MicroBenchmarkResult intersectionTest(const surface_t& surface, double phi,
                                      double theta) {
  // Shoot at it
  double cosPhi = std::cos(phi);
  double sinPhi = std::sin(phi);
  double cosTheta = std::cos(theta);
  double sinTheta = std::sin(theta);

  Vector3 direction(cosPhi * sinTheta, sinPhi * sinTheta, cosTheta);

  return Acts::Test::microBenchmark(
      [&] {
        return surface.intersect(tgContext, origin, direction, boundaryCheck);
      },
      nrepts);
}

BOOST_DATA_TEST_CASE(
    benchmark_surface_intersections,
    bdata::random(
        (bdata::seed = 21,
         bdata::distribution = std::uniform_real_distribution<>(-M_PI, M_PI))) ^
        bdata::random((bdata::seed = 22,
                       bdata::distribution =
                           std::uniform_real_distribution<>(-0.3, 0.3))) ^
        bdata::xrange(ntests),
    phi, theta, index) {
  (void)index;

  std::cout << std::endl
            << "Benchmarking theta=" << theta << ", phi=" << phi << "..."
            << std::endl;
  if (testPlane) {
    std::cout << "- Plane: "
              << intersectionTest<PlaneSurface>(*aPlane, phi, theta)
              << std::endl;
  }
  if (testDisc) {
    std::cout << "- Disc: " << intersectionTest<DiscSurface>(*aDisc, phi, theta)
              << std::endl;
  }
  if (testCylinder) {
    std::cout << "- Cylinder: "
              << intersectionTest<CylinderSurface>(*aCylinder, phi, theta)
              << std::endl;
  }
  if (testStraw) {
    std::cout << "- Straw: "
              << intersectionTest<StrawSurface>(*aStraw, phi, theta + M_PI)
              << std::endl;
  }
}

}  // namespace Test
}  // namespace Acts

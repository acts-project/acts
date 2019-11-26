// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#define BOOST_TEST_MODULE Surface Intersection Benchmark

#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
// leave blank line

#include <boost/test/data/test_case.hpp>
// leave blank line

#include <boost/test/output_test_stream.hpp>
// leave blank line

#include <cmath>
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/Surfaces/RectangleBounds.hpp"
#include "Acts/Utilities/Units.hpp"

namespace bdata = boost::unit_test::data;
namespace tt = boost::test_tools;
using namespace Acts::UnitLiterals;

namespace Acts {
namespace Test {

// Some randomness & number crunching
unsigned int ntests = 100;
unsigned int nrepts = 1000000;

// Create a test context
GeometryContext tgContext = GeometryContext();

// Create a test plane in 10 m distance
auto at = Translation3D(0., 0., 10_m);
auto rb = std::make_shared<RectangleBounds>(1_m, 1_m);
auto aPlane = Surface::makeShared<PlaneSurface>(
    std::make_shared<Transform3D>(at), std::move(rb));

Vector3D origin(0., 0., 0.);

// This test case checks that no segmentation fault appears
// - this tests the collection of surfaces
BOOST_DATA_TEST_CASE(
    test_plane_intersection,
    bdata::random(
        (bdata::seed = 21,
         bdata::distribution = std::uniform_real_distribution<>(-M_PI, M_PI))) ^
        bdata::random((bdata::seed = 22,
                       bdata::distribution =
                           std::uniform_real_distribution<>(-0.3, 0.3))) ^
        bdata::xrange(ntests),
    phi, theta, index) {
  // Shoot at it
  double cosPhi = std::cos(phi);
  double sinPhi = std::sin(phi);
  double cosTheta = std::cos(theta);
  double sinTheta = std::sin(theta);

  Vector3D direction(cosPhi * sinTheta, sinPhi * sinTheta, cosTheta);

  for (unsigned int ir = 0; ir < nrepts; ++ir) {
    auto intersect = aPlane->intersect(tgContext, origin, direction, true);

    (void)intersect;
  }
}

}  // namespace Test
}  // namespace Acts

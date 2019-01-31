// This file is part of the Acts project.
//
// Copyright (C) 2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///  Boost include(s)
#define BOOST_TEST_MODULE RelativePathCorrector Tests

#include <boost/test/included/unit_test.hpp>
// leave blank line

#include <boost/test/data/test_case.hpp>
// leave blank line

#include <boost/test/output_test_stream.hpp>
// leave blank line

#include "Acts/Propagator/detail/RelativePathCorrector.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Utilities/Helpers.hpp"

namespace bdata = boost::unit_test::data;
namespace tt    = boost::test_tools;

namespace Acts {

namespace Test {

  // This tests the implementation of the Step corrector
  BOOST_AUTO_TEST_CASE(RelativePathCorrectorTests)
  {
    // construct the surface
    Vector3D pcenter(-4., 4., 0.);
    Vector3D pnormal = Vector3D(-1., 4., 0.).normalized();
    auto     plane   = Surface::makeShared<PlaneSurface>(pcenter, pnormal);

    // The current position - origin
    Vector3D position(0., 0., 0.);
    Vector3D originalDirection = Vector3D(4., 7., 0.).normalized();
    Vector3D direction         = originalDirection;

    // The two tracks
    // (n)
    Vector3D p0n(-7., 7., 0.);
    Vector3D d0n   = Vector3D(8., 4., 0.).normalized();
    double   pathn = 10.5;
    // (p)
    Vector3D p0p(-1.5, -12., 0.);
    Vector3D d0p   = Vector3D(-2., 9., 0.).normalized();
    double   pathp = 13.;

    auto intersect = [&position, &direction, &plane]() -> Intersection {
      auto pi
          = plane->intersectionEstimate(position, direction, forward, false);
      std::cout << "Interseciton it at " << toString(pi.position) << std::endl;
      return pi;
    };

    // original intersection
    auto ointersect = intersect();

    // Step corrector for the (n) track
    detail::RelativePathCorrector corrFncn(p0n, d0n, pathn);
    // -> correct & intersect
    corrFncn(position, direction, ointersect.pathLength);
    auto nintersect = intersect();

    // re-assign the original direction
    direction = originalDirection;

    // Step corrector for the (p) track
    detail::RelativePathCorrector corrFncp(p0p, d0p, pathp);
    // -> correct & intersect
    corrFncp(position, direction, ointersect.pathLength);
    auto pintersect = intersect();

    BOOST_CHECK_LT(nintersect.pathLength, pintersect.pathLength);
  }

}  // namespace Test
}  // namespace Acts

// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/test/data/test_case.hpp>
#include <boost/test/tools/output_test_stream.hpp>
#include <boost/test/unit_test.hpp>

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Experimental/DetectorVolume.hpp"
#include "Acts/Experimental/SurfaceLinks.hpp"
#include "Acts/Geometry/CylinderVolumeBounds.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Tests/CommonHelpers/FloatComparisons.hpp"

#include <memory>
#include <vector>

#include <boost/pointer_cast.hpp>

#include "GeometryHelper.hpp"

namespace Acts {

namespace Test {

GeometryContext gctx = GeometryContext();

BOOST_AUTO_TEST_SUITE(Experimental)

BOOST_AUTO_TEST_CASE(DetectorVolumeAllSurfaces) {
  // Create a single volume with surfaces
  auto singleVolume = createBarrelVolume();
  AllSurfaces all{25};

  Vector3 atBoundary = Vector3(20., 0., 0.);
  Vector3 direction = Vector3(1., 1., 1.).normalized();

  // Ask for the surfaces that can be intersected
  auto surfaceCandidatesFromBoundary =
      all(gctx, *singleVolume, atBoundary, direction, true,
          std::numeric_limits<ActsScalar>::max(), false);

  BOOST_TEST(surfaceCandidatesFromBoundary.size(), 4u);

  // Start from the first one intersection
  auto atFirst = surfaceCandidatesFromBoundary[0].intersection.position;
  auto surfaceCandidatesFromFirst =
      all(gctx, *singleVolume, atFirst, direction, true,
          std::numeric_limits<ActsScalar>::max(), false);

  BOOST_TEST(surfaceCandidatesFromFirst.size(), 4u);

  // Start between first and boundary and go backward
  // - that will shoot to the other side
  Vector3 intermediate = 0.5 * (atBoundary + atFirst);

  auto surfaceCandidatesFromMiddle =
      all(gctx, *singleVolume, intermediate, -1 * direction, true,
          std::numeric_limits<ActsScalar>::max(), false);

  BOOST_TEST(surfaceCandidatesFromMiddle.size(), 4u);

  // Do the same, but avoid punch through
  auto surfaceCandidatesFromMiddleRestricted =
      all(gctx, *singleVolume, intermediate, -1 * direction, true, 10., false);

  BOOST_TEST(surfaceCandidatesFromMiddleRestricted.size(), 0u);
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace Test
}  // namespace Acts

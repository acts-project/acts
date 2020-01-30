// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Geometry/DoubleTrapezoidVolumeBounds.hpp"
#include "Acts/Surfaces/PlanarBounds.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Tests/CommonHelpers/FloatComparisons.hpp"
#include "Acts/Utilities/BoundingBox.hpp"
#include "Acts/Utilities/Definitions.hpp"

namespace tt = boost::test_tools;

namespace Acts {

namespace Test {
BOOST_AUTO_TEST_SUITE(Volumes)

BOOST_AUTO_TEST_CASE(bounding_box_creation) {
  float tol = 1e-4;

  DoubleTrapezoidVolumeBounds dtvb(5, 10, 8, 4, 7, 3);

  auto bb = dtvb.boundingBox();
  CHECK_CLOSE_ABS(bb.max(), Vector3D(10, 14, 3), tol);
  CHECK_CLOSE_ABS(bb.min(), Vector3D(-10, -8, -3), tol);

  Transform3D trf;

  trf = Translation3D(Vector3D(0, 30, 20));

  bb = dtvb.boundingBox(&trf);
  CHECK_CLOSE_ABS(bb.max(), Vector3D(10, 44, 23), tol);
  CHECK_CLOSE_ABS(bb.min(), Vector3D(-10, 22, 17), tol);

  trf = AngleAxis3D(M_PI / 2., Vector3D(-2, 4, 5).normalized());

  bb = dtvb.boundingBox(&trf);
  CHECK_CLOSE_ABS(bb.max(), Vector3D(8.9517, 11.7462, 10.263), tol);
  CHECK_CLOSE_ABS(bb.min(), Vector3D(-14.7572, -7.9101, -9.85174), tol);
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace Test

}  // namespace Acts

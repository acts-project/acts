// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#define BOOST_TEST_MODULE CutoutCylinderVolumeBounds Test
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <fstream>
#include <iostream>
#include <memory>

#include "Acts/Geometry/CutoutCylinderVolumeBounds.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Tests/CommonHelpers/FloatComparisons.hpp"
#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Utilities/PlyHelper.hpp"

namespace Acts {
namespace Test {

BOOST_AUTO_TEST_SUITE(Volumes)

BOOST_AUTO_TEST_CASE(construction_test) {
  CutoutCylinderVolumeBounds ccvb(5, 10, 15, 30, 25);
  ccvb.toStream(std::cout);
}

BOOST_AUTO_TEST_CASE(decomposeToSurfaces_test) {
  CutoutCylinderVolumeBounds ccvb(5, 10, 15, 30, 25);
  PlyHelper<double> ply;
  ccvb.draw(ply);

  std::ofstream os("ccvb.ply");
  os << ply;
}

BOOST_AUTO_TEST_CASE(inside_test) {
  CutoutCylinderVolumeBounds ccvb(5, 10, 15, 30, 25);

  BOOST_CHECK(!ccvb.inside({0, 0, 0}));
  BOOST_CHECK(!ccvb.inside({0, 3, 0}));
  BOOST_CHECK(!ccvb.inside({3, 0, 0}));
  BOOST_CHECK(!ccvb.inside({0, 7, 0}));
  BOOST_CHECK(!ccvb.inside({7, 0, 0}));
  BOOST_CHECK(ccvb.inside({0, 13, 0}));
  BOOST_CHECK(ccvb.inside({13, 0, 0}));
  BOOST_CHECK(!ccvb.inside({0, 17, 0}));
  BOOST_CHECK(!ccvb.inside({17, 0, 0}));

  // outside in z
  BOOST_CHECK(!ccvb.inside({0, 0, 35}));
  BOOST_CHECK(!ccvb.inside({0, 0, -35}));
  BOOST_CHECK(!ccvb.inside({0, 3, 35}));
  BOOST_CHECK(!ccvb.inside({0, 3, -35}));
  BOOST_CHECK(!ccvb.inside({3, 0, 35}));
  BOOST_CHECK(!ccvb.inside({3, 0, -35}));
  BOOST_CHECK(!ccvb.inside({0, 10, 35}));
  BOOST_CHECK(!ccvb.inside({0, 10, -35}));
  BOOST_CHECK(!ccvb.inside({10, 0, 35}));
  BOOST_CHECK(!ccvb.inside({10, 0, -35}));
  BOOST_CHECK(!ccvb.inside({0, 20, 35}));
  BOOST_CHECK(!ccvb.inside({0, 20, -35}));
  BOOST_CHECK(!ccvb.inside({20, 0, 35}));
  BOOST_CHECK(!ccvb.inside({20, 0, -35}));

  // in the choke point in z
  BOOST_CHECK(!ccvb.inside({0, 0, 27}));
  BOOST_CHECK(!ccvb.inside({0, 0, -27}));
  BOOST_CHECK(!ccvb.inside({0, 3, 27}));
  BOOST_CHECK(!ccvb.inside({0, 3, -27}));
  BOOST_CHECK(!ccvb.inside({3, 0, 27}));
  BOOST_CHECK(!ccvb.inside({3, 0, -27}));
  BOOST_CHECK(ccvb.inside({0, 7, 27}));
  BOOST_CHECK(ccvb.inside({0, 7, -27}));
  BOOST_CHECK(ccvb.inside({7, 0, 27}));
  BOOST_CHECK(ccvb.inside({7, 0, -27}));
  BOOST_CHECK(ccvb.inside({0, 13, 27}));
  BOOST_CHECK(ccvb.inside({0, 13, -27}));
  BOOST_CHECK(ccvb.inside({13, 0, 27}));
  BOOST_CHECK(ccvb.inside({13, 0, -27}));
  BOOST_CHECK(!ccvb.inside({0, 17, 27}));
  BOOST_CHECK(!ccvb.inside({0, 17, -27}));
  BOOST_CHECK(!ccvb.inside({17, 0, 27}));
  BOOST_CHECK(!ccvb.inside({17, 0, -27}));

  // right inside the choke point in z
  BOOST_CHECK(!ccvb.inside({0, 0, 23}));
  BOOST_CHECK(!ccvb.inside({0, 0, -23}));
  BOOST_CHECK(!ccvb.inside({0, 3, 23}));
  BOOST_CHECK(!ccvb.inside({0, 3, -23}));
  BOOST_CHECK(!ccvb.inside({3, 0, 23}));
  BOOST_CHECK(!ccvb.inside({3, 0, -23}));
  BOOST_CHECK(!ccvb.inside({0, 7, 23}));
  BOOST_CHECK(!ccvb.inside({0, 7, -23}));
  BOOST_CHECK(!ccvb.inside({7, 0, 23}));
  BOOST_CHECK(!ccvb.inside({7, 0, -23}));
  BOOST_CHECK(ccvb.inside({0, 13, 23}));
  BOOST_CHECK(ccvb.inside({0, 13, -23}));
  BOOST_CHECK(ccvb.inside({13, 0, 23}));
  BOOST_CHECK(ccvb.inside({13, 0, -23}));
  BOOST_CHECK(!ccvb.inside({0, 17, 23}));
  BOOST_CHECK(!ccvb.inside({0, 17, -23}));
  BOOST_CHECK(!ccvb.inside({17, 0, 23}));
  BOOST_CHECK(!ccvb.inside({17, 0, -23}));
}

BOOST_AUTO_TEST_CASE(boundingbox_test) {
  CutoutCylinderVolumeBounds ccvb(5, 10, 15, 30, 25);
  auto box = ccvb.boundingBox();
  CHECK_CLOSE_ABS(box.min(), Vector3D(-15, -15, -30), 1e-6);
  CHECK_CLOSE_ABS(box.max(), Vector3D(15, 15, 30), 1e-6);
}

BOOST_AUTO_TEST_SUITE_END()
}  // namespace Test
}  // namespace Acts

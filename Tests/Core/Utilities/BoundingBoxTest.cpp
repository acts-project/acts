// This file is part of the Acts project.
//
// Copyright (C) 2018-2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#define BOOST_TEST_MODULE BoundingBoxTests
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <fstream>
#include <iostream>
#include <map>
#include <memory>
#include <random>
#include <set>

#include "Acts/Geometry/AbstractVolume.hpp"
#include "Acts/Geometry/GenericCuboidVolumeBounds.hpp"
#include "Acts/Tests/CommonHelpers/FloatComparisons.hpp"
#include "Acts/Utilities/BoundingBox.hpp"
#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Utilities/PlyHelper.hpp"

namespace Acts {
namespace Test {

struct Object {};

using ObjectBBox = Acts::AABB3F<Object>;

BOOST_AUTO_TEST_CASE(box_construction) {
  BOOST_TEST_CONTEXT("2D") {
    Object o;
    using Box = Acts::AxisAlignedBoundingBox<Object, float, 2>;
    Box bb(&o, {-1, -1}, {2, 2});

    typename Box::transform_type rot;
    rot = Eigen::Rotation2D<float>(M_PI / 7.);
    Box bb_rot = bb.transformed(rot);

    CHECK_CLOSE_ABS(bb_rot.min(), Vector2F(-1.76874, -1.33485), 1e-4);
    CHECK_CLOSE_ABS(bb_rot.max(), Vector2F(2.23582, 2.66971), 1e-4);
  }

  BOOST_TEST_CONTEXT("3D") {
    Object o;
    using Box = Acts::AxisAlignedBoundingBox<Object, float, 3>;
    Box bb(&o, {-1, -1, -1}, {2, 2, 2});

    typename Box::transform_type rot;
    rot = AngleAxis3F(M_PI / 3., Vector3F::UnitZ());
    Box bb_rot = bb.transformed(rot);

    CHECK_CLOSE_ABS(bb_rot.min(), Vector3F(-2.23205, -1.36603, -1), 1e-4);
    CHECK_CLOSE_ABS(bb_rot.max(), Vector3F(1.86603, 2.73205, 2), 1e-4);

    rot *= AngleAxis3F(M_PI / 5., Vector3F(1, 1, 0).normalized());
    Box bb_rot2 = bb.transformed(rot);

    CHECK_CLOSE_ABS(bb_rot2.min(), Vector3F(-2.40848, -1.51816, -2.0559), 1e-4);
    CHECK_CLOSE_ABS(bb_rot2.max(), Vector3F(2.61021, 3.03631, 2.86491), 1e-4);
  }
}

BOOST_AUTO_TEST_CASE(intersect_points) {
  using vertex_type = ObjectBBox::vertex_type;

  Object o;
  ObjectBBox bb(&o, {0, 0, 0}, {1, 1, 1});
  vertex_type p;

  p = {0.5, 0.5, 0.5};
  BOOST_TEST(bb.intersect(p));
  p = {0.25, 0.25, 0.25};
  BOOST_TEST(bb.intersect(p));
  p = {0.75, 0.75, 0.75};
  BOOST_TEST(bb.intersect(p));

  // lower bound is inclusive
  p = {0, 0, 0};
  BOOST_TEST(bb.intersect(p));
  // upper bound is exclusive
  p = {1.0, 1.0, 1.0};
  BOOST_TEST(!bb.intersect(p));

  // some outsides
  p = {2, 0, 0};
  BOOST_TEST(!bb.intersect(p));
  p = {0, 2, 0};
  BOOST_TEST(!bb.intersect(p));
  p = {0, 0, 2};
  BOOST_TEST(!bb.intersect(p));
  p = {2, 2, 0};
  BOOST_TEST(!bb.intersect(p));
  p = {2, 0, 2};
  BOOST_TEST(!bb.intersect(p));
  p = {2, 2, 2};
  BOOST_TEST(!bb.intersect(p));

  p = {-1, 0, 0};
  BOOST_TEST(!bb.intersect(p));
  p = {0, -1, 0};
  BOOST_TEST(!bb.intersect(p));
  p = {0, 0, -1};
  BOOST_TEST(!bb.intersect(p));
  p = {-1, -1, 0};
  BOOST_TEST(!bb.intersect(p));
  p = {-1, 0, -1};
  BOOST_TEST(!bb.intersect(p));
  p = {-1, -1, -1};
  BOOST_TEST(!bb.intersect(p));
}

}  // namespace Test
}  // namespace Acts

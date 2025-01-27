// This file is part of the Acts project.
//
// Copyright (C) 2017-2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/test/data/test_case.hpp>
#include <boost/test/tools/output_test_stream.hpp>
#include <boost/test/unit_test.hpp>

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Surfaces/detail/VerticesHelper.hpp"

#include <algorithm>
#include <vector>

#include <Eigen/Geometry>

namespace Acts {
namespace detail {
namespace Test {
BOOST_AUTO_TEST_SUITE(Surfaces)

BOOST_AUTO_TEST_CASE(VerticesHelperOnHyperPlane) {
  // Create the transform
  Transform3 transform(AngleAxis3(0.234, Vector3(0., 1., 0.)) *
                       AngleAxis3(-0.734, Vector3(1., 1., 1.).normalized()) *
                       Translation3(Vector3(-1., 2., 3.)));

  auto trfSpace = [](std::vector<Vector3>& vtxs,
                     const Transform3& trf) -> void {
    std::transform(vtxs.begin(), vtxs.end(), vtxs.begin(),
                   [&](auto& v) { return (trf * v); });
  };

  // x-y plane test
  std::vector<Vector3> xyplane = {Vector3(1., 3., 0.), Vector3(-2., 1., 0.),
                                  Vector3(5., 8., 0.), Vector3(-9., -9., 0.),
                                  Vector3(5., 0., 0.), Vector3(3., 1., 0.)};

  trfSpace(xyplane, transform);

  // All on a hyper plane
  BOOST_CHECK(VerticesHelper::onHyperPlane(xyplane));
  // One outside the s_onSurfaceTolerance
  xyplane.push_back(transform * Vector3(3., -4., 0.05));
  BOOST_CHECK(!VerticesHelper::onHyperPlane(xyplane));
  // But inside extended tolerance
  BOOST_CHECK(VerticesHelper::onHyperPlane(xyplane, 0.6));
}

BOOST_AUTO_TEST_SUITE_END()
}  // namespace Test
}  // namespace detail
}  // namespace Acts

// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Geometry/CuboidVolumeBounds.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Tests/CommonHelpers/FloatComparisons.hpp"
#include "Acts/Utilities/Definitions.hpp"

namespace Acts {
namespace Test {

GeometryContext gctx = GeometryContext();

double hx{10.}, hy{20.}, hz{30.};

BOOST_AUTO_TEST_SUITE(Geometry)

BOOST_AUTO_TEST_CASE(CuboidVolumeConstruction) {
  // Test Construction
  CuboidVolumeBounds box(hx, hy, hz);

  // Test copy construction
  CuboidVolumeBounds copied(box);
  BOOST_CHECK_EQUAL(box, copied);

  // Test assigned
  CuboidVolumeBounds assigned = box;
  BOOST_CHECK_EQUAL(box, assigned);
}

BOOST_AUTO_TEST_CASE(CuboidVolumeRecreation) {
  CuboidVolumeBounds original(hx, hy, hz);
  auto valvector = original.values();
  std::array<double, CuboidVolumeBounds::eSize> values;
  std::copy_n(valvector.begin(), CuboidVolumeBounds::eSize, values.begin());
  CuboidVolumeBounds recreated(values);
  BOOST_CHECK_EQUAL(original, recreated);
}

BOOST_AUTO_TEST_CASE(CuboidVolumeException) {
  // Test exception negative x
  BOOST_CHECK_THROW(CuboidVolumeBounds(-hx, hy, hz), std::logic_error);
  // Test exception negative y
  BOOST_CHECK_THROW(CuboidVolumeBounds(hx, -hy, hz), std::logic_error);
  // Test exception negative z
  BOOST_CHECK_THROW(CuboidVolumeBounds(hx, hy, -hz), std::logic_error);
  // Other iterations 0
  BOOST_CHECK_THROW(CuboidVolumeBounds(-hx, hy, -hz), std::logic_error);
  // Other iterations 1
  BOOST_CHECK_THROW(CuboidVolumeBounds(-hx, -hy, hz), std::logic_error);
  // Other iterations 2
  BOOST_CHECK_THROW(CuboidVolumeBounds(hx, -hy, -hz), std::logic_error);
  // Other iterations : all
  BOOST_CHECK_THROW(CuboidVolumeBounds(-hx, -hy, -hz), std::logic_error);
}

BOOST_AUTO_TEST_CASE(CuboidVolumeProperties) {
  CuboidVolumeBounds box(hx, hy, hz);
  // Test the type
  BOOST_TEST(box.type() == VolumeBounds::eCuboid);
  // Test the halflength x
  CHECK_CLOSE_ABS(box.get(CuboidVolumeBounds::eHalfLengthX), hx, s_epsilon);
  // Test the halflength y
  CHECK_CLOSE_ABS(box.get(CuboidVolumeBounds::eHalfLengthY), hy, s_epsilon);
  // Test the halflength z
  CHECK_CLOSE_ABS(box.get(CuboidVolumeBounds::eHalfLengthZ), hz, s_epsilon);
  // Test the streaming
  std::vector<double> refvalues = {hx, hy, hz};
  BOOST_TEST(box.values() == refvalues);

  // Inside position
  Vector3D inside({5., 10., 8.});
  // Outside positions  in x, y, z
  std::vector<Vector3D> outsides = {
      {20., 1., -2.}, {1., -30., 2.}, {-1., 2., 100.}};

  // Inside position
  BOOST_TEST(box.inside(inside, s_onSurfaceTolerance));

  // Outside position
  for (const auto& outside : outsides) {
    BOOST_TEST(!box.inside(outside, s_onSurfaceTolerance));
  }
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace Test
}  // namespace Acts

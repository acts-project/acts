// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Geometry/BoundarySurfaceFace.hpp"
#include "Acts/Geometry/CuboidVolumeBounds.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Tests/CommonHelpers/FloatComparisons.hpp"

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
  std::array<double, CuboidVolumeBounds::eSize> values{};
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
  BOOST_CHECK_EQUAL(box.type(), VolumeBounds::eCuboid);
  // Test the halflength x
  CHECK_CLOSE_ABS(box.get(CuboidVolumeBounds::eHalfLengthX), hx, s_epsilon);
  // Test the halflength y
  CHECK_CLOSE_ABS(box.get(CuboidVolumeBounds::eHalfLengthY), hy, s_epsilon);
  // Test the halflength z
  CHECK_CLOSE_ABS(box.get(CuboidVolumeBounds::eHalfLengthZ), hz, s_epsilon);
  // Test the streaming
  std::vector<double> actvalues = box.values();
  std::vector<double> refvalues = {hx, hy, hz};
  BOOST_CHECK_EQUAL_COLLECTIONS(actvalues.begin(), actvalues.end(),
                                refvalues.begin(), refvalues.end());

  // Inside position
  Vector3 inside({5., 10., 8.});
  // Outside positions  in x, y, z
  std::vector<Vector3> outsides = {
      {20., 1., -2.}, {1., -30., 2.}, {-1., 2., 100.}};

  // Inside position
  BOOST_CHECK(box.inside(inside, s_onSurfaceTolerance));

  // Outside position
  for (const auto& outside : outsides) {
    BOOST_CHECK(!box.inside(outside, s_onSurfaceTolerance));
  }
}

BOOST_AUTO_TEST_CASE(CuboidVolumeBoundarySurfaces) {
  CuboidVolumeBounds box(5, 8, 7);
  auto cvbOrientedSurfaces = box.orientedSurfaces(Transform3::Identity());

  BOOST_CHECK_EQUAL(cvbOrientedSurfaces.size(), 6);

  auto geoCtx = GeometryContext();

  for (auto& os : cvbOrientedSurfaces) {
    auto osCenter = os.first->center(geoCtx);
    auto osNormal = os.first->normal(geoCtx);
    double nDir = (double)os.second;
    // Check if you step inside the volume with the oriented normal
    auto insideBox = osCenter + nDir * osNormal;
    auto outsideBox = osCenter - nDir * osNormal;
    BOOST_CHECK(box.inside(insideBox));
    BOOST_CHECK(!box.inside(outsideBox));
  }

  Vector3 xaxis(1., 0., 0.);
  Vector3 yaxis(0., 1., 0.);
  Vector3 zaxis(0., 0., 1.);

  // Test the orientation of the boundary surfaces
  auto nFaceXY =
      cvbOrientedSurfaces[negativeFaceXY].first->transform(geoCtx).rotation();
  BOOST_CHECK(nFaceXY.col(0).isApprox(xaxis));
  BOOST_CHECK(nFaceXY.col(1).isApprox(yaxis));
  BOOST_CHECK(nFaceXY.col(2).isApprox(zaxis));

  auto pFaceXY =
      cvbOrientedSurfaces[positiveFaceXY].first->transform(geoCtx).rotation();
  BOOST_CHECK(pFaceXY.col(0).isApprox(xaxis));
  BOOST_CHECK(pFaceXY.col(1).isApprox(yaxis));
  BOOST_CHECK(pFaceXY.col(2).isApprox(zaxis));

  auto nFaceYZ =
      cvbOrientedSurfaces[negativeFaceYZ].first->transform(geoCtx).rotation();
  BOOST_CHECK(nFaceYZ.col(0).isApprox(yaxis));
  BOOST_CHECK(nFaceYZ.col(1).isApprox(zaxis));
  BOOST_CHECK(nFaceYZ.col(2).isApprox(xaxis));

  auto pFaceYZ =
      cvbOrientedSurfaces[positiveFaceYZ].first->transform(geoCtx).rotation();
  BOOST_CHECK(pFaceYZ.col(0).isApprox(yaxis));
  BOOST_CHECK(pFaceYZ.col(1).isApprox(zaxis));
  BOOST_CHECK(pFaceYZ.col(2).isApprox(xaxis));

  auto nFaceZX =
      cvbOrientedSurfaces[negativeFaceZX].first->transform(geoCtx).rotation();
  BOOST_CHECK(nFaceZX.col(0).isApprox(zaxis));
  BOOST_CHECK(nFaceZX.col(1).isApprox(xaxis));
  BOOST_CHECK(nFaceZX.col(2).isApprox(yaxis));

  auto pFaceZX =
      cvbOrientedSurfaces[positiveFaceZX].first->transform(geoCtx).rotation();
  BOOST_CHECK(pFaceZX.col(0).isApprox(zaxis));
  BOOST_CHECK(pFaceZX.col(1).isApprox(xaxis));
  BOOST_CHECK(pFaceZX.col(2).isApprox(yaxis));
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace Test
}  // namespace Acts

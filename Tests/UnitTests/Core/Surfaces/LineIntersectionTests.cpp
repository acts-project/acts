// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/tools/output_test_stream.hpp>
#include <boost/test/unit_test.hpp>

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/Units.hpp"
#include "Acts/Surfaces/detail/LineHelper.hpp"
#include "Acts/Surfaces/detail/PlanarHelper.hpp"
#include "Acts/Utilities/UnitVectors.hpp"

#include <format>

using namespace Acts;

namespace ActsTests {

BOOST_AUTO_TEST_SUITE(SurfacesSuite)

BOOST_AUTO_TEST_CASE(checkPlaneIsection) {
  using namespace PlanarHelper;

  const Vector3 pointInPlane{561., 0., 0.};
  const Vector3 planeNormal{Vector3::UnitX()};

  const Vector3 extPoint = pointInPlane + 50. * planeNormal;
  const Vector3 lineDir = Vector3{1., 5., 0.}.normalized();

  const Vector3 crossPoint = intersectPlane(extPoint, lineDir, planeNormal,
                                            planeNormal.dot(pointInPlane))
                                 .position();

  BOOST_CHECK_LE(std::abs(pointInPlane.x() - crossPoint.x()),
                 std::numeric_limits<float>::epsilon());
  BOOST_CHECK_LE(std::abs(pointInPlane.y() - crossPoint.y() - 250.),
                 std::numeric_limits<float>::epsilon());
}

BOOST_AUTO_TEST_CASE(lineDistance) {
  using namespace detail::LineHelper;
  using namespace UnitLiterals;
  const Vector3 posA{100, 0., 0.};
  const Vector3 posB{100, 50, 0.};
  constexpr double tolerance = 1.e-12;

  const double parallelTest =
      signedDistance(posA, Vector3::UnitZ(), posB, Vector3::UnitZ());
  BOOST_CHECK_LE(std::abs(parallelTest - 50), tolerance);
  std::cout << __FILE__ << ":" << __LINE__
            << " - Passed test between two parallel crossing lines: "
            << parallelTest << std::endl;

  const double orthogonalTest =
      std::abs(signedDistance(posA, Vector3::UnitX(), posB, Vector3::UnitZ()));
  BOOST_CHECK_LE(std::abs(orthogonalTest - 50), tolerance);
  std::cout << __FILE__ << ":" << __LINE__
            << " - Passed test between two orthogonal lines: " << orthogonalTest
            << std::endl;

  constexpr double angleXY = 36. * 1_degree;
  const Vector3 dirInXY{makeDirectionFromPhiTheta(0., angleXY)};
  const double planeXYTest = signedDistance(posA, dirInXY, posB, dirInXY);
  BOOST_CHECK_LE(std::abs(planeXYTest - 50.), tolerance);
  const double planeXYTest1 = signedDistance(posA, dirInXY, posB, dirInXY);
  BOOST_CHECK_LE(std::abs(std::abs(planeXYTest1) - 50.), tolerance);
  std::cout << __FILE__ << ":" << __LINE__
            << " - Passed test between two parallel lines in the x-y plane: "
            << planeXYTest << std::endl;
  /// Generate a plumb-line
  const Vector3 plumbDir{-dirInXY.z(), 0., dirInXY.x()};
  /// Another vector that's perpendicular to the plumb but not parallel to the
  /// original line
  const Vector3 extDir =
      (5. * dirInXY + 21. * plumbDir.cross(dirInXY)).normalized();

  BOOST_CHECK_LE(extDir.dot(plumbDir), tolerance);
  BOOST_CHECK_LE(plumbDir.dot(dirInXY), tolerance);
  /// Create a random external point that's used as reference for the second
  /// test line
  const Vector3 extPoint =
      posA + 50. * dirInXY + 50. * plumbDir + 400. * extDir;
  /// Recuperate the external point that's closest to the primary line
  const Vector3 closePointExt =
      lineIntersect(posA, dirInXY, extPoint, extDir).position();
  BOOST_CHECK_LE((posA + 50. * dirInXY + 50. * plumbDir - closePointExt).norm(),
                 tolerance);
  /// Let's get the second closest point
  const Vector3 closePointXY =
      lineIntersect(extPoint, extDir, posA, dirInXY).position();
  BOOST_CHECK_LE((posA + 50. * dirInXY - closePointXY).norm(), tolerance);
  /// Calculate the distance between the two close points
  BOOST_CHECK_LE(std::abs((closePointXY - closePointExt).norm() - 50.),
                 tolerance);

  const double extLineDist = signedDistance(posA, dirInXY, extPoint, extDir);
  BOOST_CHECK_LE(std::abs(extLineDist - 50.), tolerance);
  const double extLineDist1 = signedDistance(posA, dirInXY, extPoint, extDir);
  BOOST_CHECK_LE(std::abs(extLineDist1 - 50.), tolerance);
  /// Finally the case when both lines are crossing
  const double crossLines =
      signedDistance(extPoint + 525. * dirInXY, dirInXY, extPoint, extDir);
  // Use a larger tolerance here.  This is the result of the sqrt of the
  // difference of numbers O(1e-5), so even with double precision, we can't
  // expect much better than this.
  BOOST_CHECK_LE(std::abs(crossLines), 1e-5);
}
BOOST_AUTO_TEST_SUITE_END()

}  // namespace ActsTests

// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/Units.hpp"
#include "Acts/Geometry/CuboidVolumeBounds.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/Volume.hpp"
#include "Acts/Utilities/BinningType.hpp"
#include "ActsTests/CommonHelpers/FloatComparisons.hpp"

#include <cmath>
#include <limits>
#include <memory>
#include <utility>

using namespace Acts;

namespace ActsTests {

BOOST_AUTO_TEST_SUITE(GeometrySuite)

BOOST_AUTO_TEST_CASE(VolumeTest) {
  using namespace UnitLiterals;
  double eps = std::numeric_limits<double>::epsilon();

  const auto gctx = GeometryContext::dangerouslyDefaultConstruct();

  // Build a translation
  Vector3 translation{1_mm, 2_mm, 3_mm};

  // Build a translation
  SquareMatrix3 rotation = RotationMatrix3::Identity();
  double rotationAngle = 60_degree;
  Vector3 xPos(std::cos(rotationAngle), 0., std::sin(rotationAngle));
  Vector3 yPos(0., 1., 0.);
  Vector3 zPos(-std::sin(rotationAngle), 0., std::cos(rotationAngle));
  rotation.col(0) = xPos;
  rotation.col(1) = yPos;
  rotation.col(2) = zPos;

  // Build a transform
  Transform3 transform(Transform3::Identity() * rotation);
  transform.translation() = translation;
  // Build the bounds
  CuboidVolumeBounds bounds(4_mm, 5_mm, 6_mm);

  // Build and test the volume
  Volume volume(transform, std::make_shared<CuboidVolumeBounds>(bounds));
  BOOST_CHECK_EQUAL(volume.localToGlobalTransform(gctx).matrix(),
                    transform.matrix());
  CHECK_CLOSE_ABS(volume.globalToLocalTransform(gctx).matrix(),
                  transform.inverse().matrix(), eps);
  BOOST_CHECK_EQUAL(volume.center(gctx), translation);
  auto vBounds = static_cast<const decltype(bounds)*>(&volume.volumeBounds());
  BOOST_CHECK_EQUAL(*vBounds, bounds);

  // Build and test a shifted volume
  Transform3 shift(Transform3::Identity());
  Vector3 shiftTranslation{-4_mm, -5_mm, -6_mm};
  shift.translation() = shiftTranslation;
  Volume volumeShift = volume.shifted(gctx, shift);
  BOOST_CHECK_EQUAL(
      volumeShift.center(gctx),
      (shift * volume.localToGlobalTransform(gctx)).translation());
  BOOST_CHECK_EQUAL(volumeShift.localToGlobalTransform(gctx).rotation(),
                    volume.localToGlobalTransform(gctx).rotation());

  // Inside/Outside check
  BOOST_CHECK(volume.inside(gctx, translation));
  BOOST_CHECK(!volume.inside(gctx, {10_mm, 2_mm, 3_mm}));
  BOOST_CHECK(volume.inside(gctx, {10_mm, 2_mm, 3_mm}, 2_mm));
  BOOST_CHECK_EQUAL(volume.referencePosition(gctx, AxisDirection::AxisX),
                    volume.center(gctx));
}

BOOST_AUTO_TEST_CASE(VolumeUpdateTest) {
  using namespace UnitLiterals;
  auto volBounds = std::make_shared<CuboidVolumeBounds>(4_mm, 5_mm, 6_mm);
  auto volBounds2 = std::make_shared<CuboidVolumeBounds>(4_mm, 5_mm, 8_mm);

  const auto gctx = GeometryContext::dangerouslyDefaultConstruct();

  Transform3 trf = Transform3::Identity();

  Volume volume(trf, volBounds);

  // Only update the bounds, keep the transform the same
  volume.update(gctx, volBounds2, std::nullopt);
  BOOST_CHECK_EQUAL(&volume.volumeBounds(), volBounds2.get());
  BOOST_CHECK_EQUAL(volume.localToGlobalTransform(gctx).matrix(), trf.matrix());

  // Update the bounds and the transform
  Transform3 trf2{Translation3{1_mm, 2_mm, 3_mm}};
  volume.update(gctx, volBounds, trf2);
  BOOST_CHECK_EQUAL(&volume.volumeBounds(), volBounds.get());
  BOOST_CHECK_EQUAL(volume.localToGlobalTransform(gctx).matrix(),
                    trf2.matrix());
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace ActsTests

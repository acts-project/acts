// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/test/tools/old/interface.hpp>
#include <boost/test/unit_test.hpp>
#include <boost/test/unit_test_suite.hpp>

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/Units.hpp"
#include "Acts/Geometry/CuboidVolumeBounds.hpp"
#include "Acts/Geometry/CylinderVolumeBounds.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/Volume.hpp"
#include "Acts/Tests/CommonHelpers/FloatComparisons.hpp"
#include "Acts/Utilities/BinningType.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <memory>

using namespace Acts::UnitLiterals;

namespace Acts::Test {

BOOST_AUTO_TEST_SUITE(Geometry);

BOOST_AUTO_TEST_CASE(VolumeTest) {
  double eps = std::numeric_limits<double>::epsilon();

  // Build a translation
  Vector3 translation{1_mm, 2_mm, 3_mm};

  // Build a translation
  ActsMatrix<3, 3> rotation = RotationMatrix3::Identity();
  double rotationAngle = 60_degree;
  Vector3 xPos(cos(rotationAngle), 0., sin(rotationAngle));
  Vector3 yPos(0., 1., 0.);
  Vector3 zPos(-sin(rotationAngle), 0., cos(rotationAngle));
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
  BOOST_CHECK_EQUAL(volume.transform().matrix(), transform.matrix());
  CHECK_CLOSE_ABS(volume.itransform().matrix(), transform.inverse().matrix(),
                  eps);
  BOOST_CHECK_EQUAL(volume.center(), translation);
  auto vBounds = static_cast<const decltype(bounds)*>(&volume.volumeBounds());
  BOOST_CHECK_EQUAL(*vBounds, bounds);

  // Build and test a shifted volume
  Transform3 shift(Transform3::Identity());
  Vector3 shiftTranslation{-4_mm, -5_mm, -6_mm};
  shift.translation() = shiftTranslation;
  Volume volumeShift(volume, shift);
  BOOST_CHECK_EQUAL(volumeShift.center(),
                    (shift * volume.transform()).translation());
  BOOST_CHECK_EQUAL(volumeShift.transform().rotation(),
                    volume.transform().rotation());

  // Inside/Outside check
  BOOST_CHECK(volume.inside(translation));
  BOOST_CHECK(!volume.inside({10_mm, 2_mm, 3_mm}));
  BOOST_CHECK(volume.inside({10_mm, 2_mm, 3_mm}, 2_mm));

  // Binning test
  GeometryContext gctx;
  BOOST_CHECK_EQUAL(volume.binningPosition(gctx, binX), volume.center());
}

BOOST_AUTO_TEST_CASE(VolumeUpdateTest) {
  using namespace Acts::UnitLiterals;
  auto volBounds = std::make_shared<CuboidVolumeBounds>(4_mm, 5_mm, 6_mm);
  auto volBounds2 = std::make_shared<CuboidVolumeBounds>(4_mm, 5_mm, 8_mm);

  Transform3 trf = Transform3::Identity();

  Volume volume(trf, volBounds);

  // Only update the bounds, keep the transform the same
  volume.update(volBounds2, std::nullopt);
  BOOST_CHECK_EQUAL(&volume.volumeBounds(), volBounds2.get());
  BOOST_CHECK_EQUAL(volume.transform().matrix(), trf.matrix());

  // Update the bounds and the transform
  Transform3 trf2{Translation3{1_mm, 2_mm, 3_mm}};
  volume.update(volBounds, trf2);
  BOOST_CHECK_EQUAL(&volume.volumeBounds(), volBounds.get());
  BOOST_CHECK_EQUAL(volume.transform().matrix(), trf2.matrix());
}

std::vector<std::shared_ptr<Volume>> joinCylinderVolumes(
    const std::vector<std::shared_ptr<Volume>>& volumes,
    BinningValue direction) {
  std::vector<CylinderVolumeBounds*> cylinderBounds;
  std::transform(
      volumes.begin(), volumes.end(), std::back_inserter(cylinderBounds),
      [](auto& vol) -> CylinderVolumeBounds* {
        return dynamic_cast<CylinderVolumeBounds*>(&vol->volumeBounds());
      });

  if (std::any_of(cylinderBounds.begin(), cylinderBounds.end(),
                  [](auto& bounds) { return bounds == nullptr; })) {
    throw std::invalid_argument("Not all volumes are cylinder volumes");
  }

  if (direction == Acts::binZ) {
    // @TODO Check if the volumes are aligned

    // Synchronize in R: collect global min and max
    ActsScalar minR =
        (*std::min_element(cylinderBounds.begin(), cylinderBounds.end(),
                           [](auto* a, auto* b) {
                             return a->get(CylinderVolumeBounds::eMinR) <
                                    b->get(CylinderVolumeBounds::eMinR);
                           }))
            ->get(CylinderVolumeBounds::eMinR);

    ActsScalar maxR =
        (*std::max_element(cylinderBounds.begin(), cylinderBounds.end(),
                           [](auto* a, auto* b) {
                             return a->get(CylinderVolumeBounds::eMaxR) <
                                    b->get(CylinderVolumeBounds::eMaxR);
                           }))
            ->get(CylinderVolumeBounds::eMaxR);

    for (auto* bounds : cylinderBounds) {
      bounds->set({
          {CylinderVolumeBounds::eMinR, minR},
          {CylinderVolumeBounds::eMaxR, maxR},
      });
    }

  } else {
    throw std::invalid_argument(binningValueNames()[direction] +
                                " is not supported ");
  }

  return volumes;
}

BOOST_AUTO_TEST_CASE(JoinCylinderVolumesAlongZ) {
  // @TODO: Rotated cylinder volumes

  ActsScalar hlZ = 400_mm;

  // Cylinder volumes which already line up, but have different1 radii
  auto bounds1 = std::make_shared<CylinderVolumeBounds>(100_mm, 400_mm, hlZ);
  auto bounds2 = std::make_shared<CylinderVolumeBounds>(200_mm, 600_mm, hlZ);
  auto bounds3 = std::make_shared<CylinderVolumeBounds>(300_mm, 500_mm, hlZ);

  auto transform = Transform3::Identity();

  transform.translation() << 0_mm, 0_mm, -hlZ;
  auto vol1 = std::make_shared<Volume>(transform, bounds1);

  transform.translation() << 0_mm, 0_mm, 0_mm;
  auto vol2 = std::make_shared<Volume>(transform, bounds2);

  transform.translation() << 0_mm, 0_mm, hlZ;
  auto vol3 = std::make_shared<Volume>(transform, bounds3);

  std::vector<std::shared_ptr<Volume>> volumes = {vol1, vol2, vol3};
  auto joined = joinCylinderVolumes(volumes, binZ);

  BOOST_CHECK_EQUAL(joined.size(), 3);

  // We get back the original volumes
  BOOST_CHECK_EQUAL_COLLECTIONS(joined.begin(), joined.end(), volumes.begin(),
                                volumes.end());

  std::vector<CylinderVolumeBounds*> origBounds = {bounds1.get(), bounds2.get(),
                                                   bounds3.get()};
  std::vector<VolumeBounds*> newBounds = {&joined[0]->volumeBounds(),
                                          &joined[1]->volumeBounds(),
                                          &joined[2]->volumeBounds()};

  // The original bounds were modified
  BOOST_CHECK_EQUAL_COLLECTIONS(origBounds.begin(), origBounds.end(),
                                newBounds.begin(), newBounds.end());

  for (const auto* bounds : origBounds) {
    BOOST_CHECK_EQUAL(bounds->get(CylinderVolumeBounds::eMinR), 100_mm);
    BOOST_CHECK_EQUAL(bounds->get(CylinderVolumeBounds::eMaxR), 600_mm);
  }
}

BOOST_AUTO_TEST_SUITE_END();

}  // namespace Acts::Test

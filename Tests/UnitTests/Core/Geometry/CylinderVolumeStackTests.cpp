// This file is part of the Acts project.
//
// Copyright (C) 2024 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/test/data/test_case.hpp>
#include <boost/test/tools/old/interface.hpp>
#include <boost/test/unit_test.hpp>
#include <boost/test/unit_test_suite.hpp>

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/Units.hpp"
#include "Acts/Geometry/CylinderVolumeBounds.hpp"
#include "Acts/Geometry/CylinderVolumeStack.hpp"
#include "Acts/Tests/CommonHelpers/FloatComparisons.hpp"

using namespace Acts::UnitLiterals;

auto logger = Acts::getDefaultLogger("UnitTests", Acts::Logging::VERBOSE);

namespace Acts::Test {

BOOST_AUTO_TEST_SUITE(Geometry);

// @TODO: Rotated cylinder volumes

BOOST_DATA_TEST_CASE(
    JoinCylinderVolumesAlongZ,
    boost::unit_test::data::xrange(-135, 180, 45) *
        boost::unit_test::data::make(0.8, 1.0, 1.2) *
        boost::unit_test::data::make(Vector3{0_mm, 0_mm, 0_mm},
                                     Vector3{20_mm, 0_mm, 0_mm},
                                     Vector3{0_mm, 20_mm, 0_mm},
                                     Vector3{20_mm, 20_mm, 0_mm},
                                     Vector3{0_mm, 0_mm, 20_mm}),
    angle, shift, offset) {
  ActsScalar hlZ = 400_mm;

  // Cylinder volumes which already line up, but have different1 radii
  auto bounds1 = std::make_shared<CylinderVolumeBounds>(100_mm, 400_mm, hlZ);
  auto bounds2 = std::make_shared<CylinderVolumeBounds>(200_mm, 600_mm, hlZ);
  auto bounds3 = std::make_shared<CylinderVolumeBounds>(300_mm, 500_mm, hlZ);

  Transform3 base =
      AngleAxis3(angle * 1_degree, Vector3::UnitX()) * Translation3(offset);

  Transform3 transform1 = base;
  transform1.translate(Vector3{0_mm, 0_mm, -2 * hlZ * shift});
  auto vol1 = std::make_shared<Volume>(transform1, bounds1);

  Transform3 transform2 = base;
  transform2.translate(Vector3{0_mm, 0_mm, 0_mm});
  auto vol2 = std::make_shared<Volume>(transform2, bounds2);

  Transform3 transform3 = base;
  transform3.translate(Vector3{0_mm, 0_mm, 2 * hlZ * shift});
  auto vol3 = std::make_shared<Volume>(transform3, bounds3);

  std::vector<std::shared_ptr<Volume>> volumes = {vol1, vol2, vol3};

  if (shift < 1.0) {
    BOOST_CHECK_THROW(
        CylinderVolumeStack(volumes, binZ,
                            CylinderVolumeStack::AttachmentStrategy::Gap,
                            *logger),
        std::invalid_argument);
    return;
  }
  auto cylStack = CylinderVolumeStack(
      volumes, binZ, CylinderVolumeStack::AttachmentStrategy::Gap, *logger);

  auto stackBounds =
      dynamic_cast<const CylinderVolumeBounds*>(&cylStack.volumeBounds());
  BOOST_REQUIRE(stackBounds != nullptr);

  BOOST_CHECK_EQUAL(stackBounds->get(CylinderVolumeBounds::eMinR), 100_mm);
  BOOST_CHECK_EQUAL(stackBounds->get(CylinderVolumeBounds::eMaxR), 600_mm);
  BOOST_CHECK_EQUAL(stackBounds->get(CylinderVolumeBounds::eHalfLengthZ),
                    hlZ + 2 * hlZ * shift);
  CHECK_CLOSE_OR_SMALL(cylStack.transform().matrix(), base.matrix(), 1e-10,
                       1e-14);

  // No expansion, original volumes did not move
  BOOST_CHECK_EQUAL(vol1->transform().matrix(), transform1.matrix());
  BOOST_CHECK_EQUAL(vol2->transform().matrix(), transform2.matrix());
  BOOST_CHECK_EQUAL(vol3->transform().matrix(), transform3.matrix());

  // All volumes (including gaps) are cylinders and have the same radial bounds
  for (const auto& volume : volumes) {
    const auto* cylinderBounds =
        dynamic_cast<const CylinderVolumeBounds*>(&volume->volumeBounds());
    BOOST_REQUIRE(cylinderBounds != nullptr);
    BOOST_CHECK_EQUAL(cylinderBounds->get(CylinderVolumeBounds::eMinR), 100_mm);
    BOOST_CHECK_EQUAL(cylinderBounds->get(CylinderVolumeBounds::eMaxR), 600_mm);
  }

  if (shift <= 1.0) {
    // No gap volumes were added
    BOOST_CHECK_EQUAL(volumes.size(), 3);
  } else {
    // Gap volumes were added
    BOOST_CHECK_EQUAL(volumes.size(), 5);
  }
}

// @TODO: Merge z asymmetric
// @TODO: Zero volumes
// @TODO: Single volumes

#if 0
BOOST_AUTO_TEST_CASE(JoinCylinderVolumesAlongZNoGapsCommonRotation) {
  ActsScalar hlZ = 400_mm;

  // Cylinder volumes which already line up, but have different1 radii
  auto bounds1 = std::make_shared<CylinderVolumeBounds>(100_mm, 400_mm, hlZ);
  auto bounds2 = std::make_shared<CylinderVolumeBounds>(200_mm, 600_mm, hlZ);
  auto bounds3 = std::make_shared<CylinderVolumeBounds>(300_mm, 500_mm, hlZ);

  auto transform = Transform3::Identity();

  auto rot = AngleAxis3(90_degree, Vector3::UnitX());

  transform.translation() << 0_mm, 0_mm, -hlZ;
  transform.prerotate(rot);
  auto vol1 = std::make_shared<Volume>(transform, bounds1);

  transform.translation() << 0_mm, 0_mm, 0_mm;
  transform.prerotate(rot);
  auto vol2 = std::make_shared<Volume>(transform, bounds2);

  transform.translation() << 0_mm, 0_mm, hlZ;
  transform.prerotate(rot);
  auto vol3 = std::make_shared<Volume>(transform, bounds3);

  // for (const auto& vol : {vol1, vol2, vol3}) {
  // std::cout << vol->transform().matrix() << std::endl;
  // }
}

BOOST_AUTO_TEST_CASE(JoinCylinderVolumesAlongZMisaligned) {
  ActsScalar hlZ = 400_mm;

  auto bounds1 = std::make_shared<CylinderVolumeBounds>(100_mm, 400_mm, hlZ);
  auto bounds2 = std::make_shared<CylinderVolumeBounds>(200_mm, 600_mm, hlZ);
  auto bounds3 = std::make_shared<CylinderVolumeBounds>(300_mm, 500_mm, hlZ);

  auto transform = Transform3::Identity();

  {
    // rotated differently

    transform.setIdentity();
    transform.translation() << 0_mm, 0_mm, -hlZ;
    auto vol1 = std::make_shared<Volume>(transform, bounds1);

    transform.setIdentity();
    transform.rotate(Eigen::AngleAxis(30_degree, Eigen::Vector3d::UnitX()));
    transform.translation() << 0_mm, 0_mm, 0_mm;
    auto vol2 = std::make_shared<Volume>(transform, bounds2);

    transform.setIdentity();
    transform.rotate(Eigen::AngleAxis(30_degree, Eigen::Vector3d::UnitY()));
    transform.translation() << 0_mm, 0_mm, hlZ;
    auto vol3 = std::make_shared<Volume>(transform, bounds3);

    BOOST_CHECK_THROW(joinCylinderVolumes({vol1, vol2, vol3}, binZ),
                      std::invalid_argument);
  }

  {
    // aligned with global z but shifted xy

    transform.setIdentity();
    transform.translation() << -10_mm, 0_mm, -hlZ;
    auto vol1 = std::make_shared<Volume>(transform, bounds1);

    transform.setIdentity();
    transform.translation() << 0_mm, 10_mm, 0_mm;
    auto vol2 = std::make_shared<Volume>(transform, bounds2);

    transform.setIdentity();
    transform.translation() << 0_mm, 0_mm, hlZ;
    auto vol3 = std::make_shared<Volume>(transform, bounds3);

    BOOST_CHECK_THROW(joinCylinderVolumes({vol1, vol2, vol3}, binZ),
                      std::invalid_argument);
  }

  {
    // aligned with global z but shifted xy

    transform.setIdentity();
    transform.translation() << -10_mm, 0_mm, -hlZ;
    auto vol1 = std::make_shared<Volume>(transform, bounds1);

    transform.setIdentity();
    transform.translation() << 0_mm, 10_mm, 0_mm;
    auto vol2 = std::make_shared<Volume>(transform, bounds2);

    transform.setIdentity();
    transform.translation() << 0_mm, 0_mm, hlZ;
    auto vol3 = std::make_shared<Volume>(transform, bounds3);

    BOOST_CHECK_THROW(joinCylinderVolumes({vol1, vol2, vol3}, binZ),
                      std::invalid_argument);
  }
}

BOOST_AUTO_TEST_CASE(JoinCylinderVolumesAlongZOverlapping) {
  ActsScalar hlZ = 400_mm;

  auto bounds1 = std::make_shared<CylinderVolumeBounds>(100_mm, 400_mm, hlZ);
  auto bounds2 = std::make_shared<CylinderVolumeBounds>(200_mm, 600_mm, hlZ);
  auto bounds3 = std::make_shared<CylinderVolumeBounds>(300_mm, 500_mm, hlZ);

  auto transform = Transform3::Identity();

  transform.translation() << 0_mm, 0_mm, -hlZ * 0.7;
  auto vol1 = std::make_shared<Volume>(transform, bounds1);

  transform.translation() << 0_mm, 0_mm, 0_mm;
  auto vol2 = std::make_shared<Volume>(transform, bounds2);

  transform.translation() << 0_mm, 0_mm, hlZ;
  auto vol3 = std::make_shared<Volume>(transform, bounds3);

  BOOST_CHECK_THROW(joinCylinderVolumes({vol1, vol2, vol3}, binZ),
                    std::invalid_argument);
}

BOOST_AUTO_TEST_CASE(JoinCylinderVolumesAlongZWithGaps) {
  ActsScalar hlZ = 400_mm;

  auto bounds1 = std::make_shared<CylinderVolumeBounds>(100_mm, 400_mm, hlZ);
  auto bounds2 = std::make_shared<CylinderVolumeBounds>(200_mm, 600_mm, hlZ);
  auto bounds3 = std::make_shared<CylinderVolumeBounds>(300_mm, 500_mm, hlZ);

  auto transform = Transform3::Identity();

  transform.translation() << 0_mm, 0_mm, -hlZ * 1.2;
  auto vol1 = std::make_shared<Volume>(transform, bounds1);

  transform.translation() << 0_mm, 0_mm, 0_mm;
  auto vol2 = std::make_shared<Volume>(transform, bounds2);

  transform.translation() << 0_mm, 0_mm, hlZ * 1.2;
  auto vol3 = std::make_shared<Volume>(transform, bounds3);

  std::vector<std::shared_ptr<Volume>> volumes = {vol1, vol2, vol3};
  auto joined = joinCylinderVolumes(volumes, binZ);

  BOOST_CHECK_EQUAL(joined.size(), 5);
  // These three are the original volumes
  BOOST_CHECK_EQUAL(joined.at(0), vol1);
  BOOST_CHECK_EQUAL(joined.at(2), vol2);
  BOOST_CHECK_EQUAL(joined.at(4), vol3);
}

#endif

BOOST_AUTO_TEST_SUITE_END();

}  // namespace Acts::Test

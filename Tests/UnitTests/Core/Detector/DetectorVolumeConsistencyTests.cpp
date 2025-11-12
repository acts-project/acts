// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Detector/DetectorVolume.hpp"
#include "Acts/Detector/PortalGenerators.hpp"
#include "Acts/Detector/detail/DetectorVolumeConsistency.hpp"
#include "Acts/Geometry/CuboidVolumeBounds.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Navigation/InternalNavigation.hpp"

using namespace Acts;

auto portalGenerator = Experimental::defaultPortalGenerator();
auto tContext = GeometryContext();

namespace ActsTests {

BOOST_AUTO_TEST_SUITE(DetectorSuite)

BOOST_AUTO_TEST_CASE(DetectorVolumeConsistencyFail) {
  // A perfect box shape
  auto box = std::make_shared<CuboidVolumeBounds>(10, 10, 10);

  // Create volume A
  auto volumeA = Experimental::DetectorVolumeFactory::construct(
      portalGenerator, tContext, "VolumeA", Transform3::Identity(), box,
      Experimental::tryAllPortals());

  // Move it into the bval direction
  auto transformB = Transform3::Identity();
  Vector3 translationB = Vector3::Zero();
  translationB[toUnderlying(AxisDirection::AxisX)] = 20;
  translationB[toUnderlying(AxisDirection::AxisY)] = 5;
  transformB.pretranslate(translationB);
  // Create volume B
  auto volumeB = Experimental::DetectorVolumeFactory::construct(
      portalGenerator, tContext, "VolumeB", transformB, box,
      Experimental::tryAllPortals());
  // Build the container
  std::vector<std::shared_ptr<Experimental::DetectorVolume>> volumes = {
      volumeA, volumeB};

  BOOST_CHECK_THROW(
      Experimental::detail::DetectorVolumeConsistency::checkCenterAlignment(
          tContext, {volumeA, volumeB}, AxisDirection::AxisX),
      std::invalid_argument);
}

BOOST_AUTO_TEST_CASE(DetectorVolumeConsistencyPass) {
  // A perfect box shape
  auto box = std::make_shared<CuboidVolumeBounds>(10, 10, 10);

  // Create volume A
  auto volumeA = Experimental::DetectorVolumeFactory::construct(
      portalGenerator, tContext, "VolumeA", Transform3::Identity(), box,
      Experimental::tryAllPortals());

  // Move it into the bval direction
  auto transformB = Transform3::Identity();
  Vector3 translationB = Vector3::Zero();
  translationB[toUnderlying(AxisDirection::AxisX)] = 20;
  transformB.pretranslate(translationB);
  // Create volume B
  auto volumeB = Experimental::DetectorVolumeFactory::construct(
      portalGenerator, tContext, "VolumeB", transformB, box,
      Experimental::tryAllPortals());
  // Build the container
  std::vector<std::shared_ptr<Experimental::DetectorVolume>> volumes = {
      volumeA, volumeB};

  BOOST_CHECK_NO_THROW(
      Experimental::detail::DetectorVolumeConsistency::checkCenterAlignment(
          tContext, {volumeA, volumeB}, AxisDirection::AxisX));
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace ActsTests

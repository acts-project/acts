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
#include "Acts/Detector/detail/CuboidalDetectorHelper.hpp"
#include "Acts/Geometry/CuboidVolumeBounds.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Navigation/InternalNavigation.hpp"
#include "Acts/Utilities/BinningData.hpp"
#include "Acts/Utilities/StringHelpers.hpp"
#include "Acts/Visualization/GeometryView3D.hpp"
#include "Acts/Visualization/ObjVisualization3D.hpp"

#include <array>
#include <memory>
#include <stdexcept>

using namespace Acts;

auto portalGenerator = Experimental::defaultPortalGenerator();
auto tContext = GeometryContext();

namespace ActsTests {

BOOST_AUTO_TEST_SUITE(DetectorSuite)

BOOST_AUTO_TEST_CASE(CubicVolumeExceptions) {
  // A perfect box shape
  auto box = std::make_shared<CuboidVolumeBounds>(10, 10, 10);

  // Create volume A
  auto volumeA = Experimental::DetectorVolumeFactory::construct(
      portalGenerator, tContext, "VolumeA", Transform3::Identity(), box,
      Experimental::tryAllPortals());

  // Create volume B
  auto transformB = Transform3::Identity();
  transformB.prerotate(AngleAxis3(0.234, Vector3::UnitZ()));

  auto volumeB = Experimental::DetectorVolumeFactory::construct(
      portalGenerator, tContext, "volumeB", transformB, box,
      Experimental::tryAllPortals());
  // Build the container
  std::vector<std::shared_ptr<Experimental::DetectorVolume>> volumes = {
      volumeA, volumeB};

  BOOST_CHECK_THROW(
      Experimental::detail::CuboidalDetectorHelper::connect(
          tContext, volumes, AxisDirection::AxisX, {}, Logging::VERBOSE),
      std::invalid_argument);
}

BOOST_AUTO_TEST_CASE(SimpleBoxConnection) {
  std::array<AxisDirection, 3> binningValues = {
      AxisDirection::AxisX, AxisDirection::AxisY, AxisDirection::AxisZ};
  for (auto bVal : binningValues) {
    // A perfect box shape
    auto box = std::make_shared<CuboidVolumeBounds>(10, 10, 10);

    // Create volume A
    auto volumeA = Experimental::DetectorVolumeFactory::construct(
        portalGenerator, tContext, "VolumeA", Transform3::Identity(), box,
        Experimental::tryAllPortals());

    // Move it into the bval direction
    auto transformB = Transform3::Identity();

    Vector3 translation = Vector3::Zero();
    translation[toUnderlying(bVal)] = 20;
    transformB.pretranslate(translation);
    // Create volume B
    auto volumeB = Experimental::DetectorVolumeFactory::construct(
        portalGenerator, tContext, "VolumeB", transformB, box,
        Experimental::tryAllPortals());
    // Build the container
    std::vector<std::shared_ptr<Experimental::DetectorVolume>> volumes = {
        volumeA, volumeB};
    auto container = Experimental::detail::CuboidalDetectorHelper::connect(
        tContext, volumes, bVal, {}, Logging::VERBOSE);

    // Check the container is constructed
    BOOST_CHECK(!container.empty());

    ObjVisualization3D obj;
    GeometryView3D::drawDetectorVolume(obj, *volumeA, tContext);
    GeometryView3D::drawDetectorVolume(obj, *volumeB, tContext);
    obj.write("ConnectectBoxesRegular_" + axisDirectionName(bVal) + ".obj");
  }
}

BOOST_AUTO_TEST_CASE(IrregularBoxConnectionInZ) {
  std::vector<AxisDirection> binningValues = {
      AxisDirection::AxisX, AxisDirection::AxisY, AxisDirection::AxisZ};

  using HlPos = std::array<double, 2u>;
  using VolHlPos = std::array<HlPos, 3u>;
  using VolSetup = std::array<VolHlPos, 3u>;

  VolHlPos cPA = {{{10., 0.}, {10., 0.}, {10., 0.}}};
  VolHlPos cPB = {{{20., 0.}, {20., 0.}, {20., 0.}}};
  VolHlPos sP = {{{10., -30.}, {30., 10.}, {90., 130.}}};

  std::array<VolSetup, 3u> volSetups = {
      {{sP, cPA, cPB}, {cPB, sP, cPA}, {cPA, cPB, sP}}};

  std::array<Transform3, 2u> transforms = {
      Transform3::Identity(),
      Transform3{Transform3::Identity()}.prerotate(
          AngleAxis3(0.34, Vector3(1., 1., 1.).normalized()))};

  // Try with arbitrary rotations
  for (auto [it, t] : enumerate(transforms)) {
    std::string trstr = it == 0 ? "" : "_rotated";
    auto rotation = t.rotation();
    // Try for all binning values
    for (auto bVal : binningValues) {
      auto [vsA, vsB, vsC] = volSetups[toUnderlying(bVal)];

      // Three box shares with different length in Z
      auto boxA =
          std::make_shared<CuboidVolumeBounds>(vsA[0][0], vsB[0][0], vsC[0][0]);
      auto boxB =
          std::make_shared<CuboidVolumeBounds>(vsA[1][0], vsB[1][0], vsC[1][0]);
      auto boxC =
          std::make_shared<CuboidVolumeBounds>(vsA[2][0], vsB[2][0], vsC[2][0]);

      auto transformA = Transform3::Identity();
      auto transformB = Transform3::Identity();
      auto transformC = Transform3::Identity();

      transformA.pretranslate(Vector3(vsA[0][1], vsB[0][1], vsC[0][1]));
      transformB.pretranslate(Vector3(vsA[1][1], vsB[1][1], vsC[1][1]));
      transformC.pretranslate(Vector3(vsA[2][1], vsB[2][1], vsC[2][1]));

      transformA.prerotate(rotation);
      transformB.prerotate(rotation);
      transformC.prerotate(rotation);

      // Create volume A, B, C
      auto volumeA = Experimental::DetectorVolumeFactory::construct(
          portalGenerator, tContext, "VolumeA", transformA, boxA,
          Experimental::tryAllPortals());
      auto volumeB = Experimental::DetectorVolumeFactory::construct(
          portalGenerator, tContext, "VolumeB", transformB, boxB,
          Experimental::tryAllPortals());
      auto volumeC = Experimental::DetectorVolumeFactory::construct(
          portalGenerator, tContext, "VolumeC", transformC, boxC,
          Experimental::tryAllPortals());

      // Build the container
      std::vector<std::shared_ptr<Experimental::DetectorVolume>> volumes = {
          volumeA, volumeB, volumeC};
      auto container = Experimental::detail::CuboidalDetectorHelper::connect(
          tContext, volumes, bVal, {}, Logging::VERBOSE);

      // Check the container is constructed
      BOOST_CHECK(!container.empty());

      ObjVisualization3D obj;
      GeometryView3D::drawDetectorVolume(obj, *volumeA, tContext);
      GeometryView3D::drawDetectorVolume(obj, *volumeB, tContext);
      GeometryView3D::drawDetectorVolume(obj, *volumeC, tContext);
      obj.write("ConnectectBoxesIrregular_" + axisDirectionName(bVal) + trstr +
                ".obj");
    }
  }
}

BOOST_AUTO_TEST_CASE(ContainerConnection) {
  // A perfect box shape
  auto box = std::make_shared<CuboidVolumeBounds>(10, 10, 10);

  // Create an AB container

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
  auto containerAB = Experimental::detail::CuboidalDetectorHelper::connect(
      tContext, volumes, AxisDirection::AxisX, {}, Logging::VERBOSE);

  // Create a CD container

  auto transformC = Transform3::Identity();
  Vector3 translationC = Vector3::Zero();
  translationC[toUnderlying(AxisDirection::AxisY)] = 20;
  transformC.pretranslate(translationC);

  auto volumeC = Experimental::DetectorVolumeFactory::construct(
      portalGenerator, tContext, "VolumeC", transformC, box,
      Experimental::tryAllPortals());

  auto transformD = Transform3::Identity();
  Vector3 translationD = Vector3::Zero();
  translationD[toUnderlying(AxisDirection::AxisX)] = 20;
  translationD[toUnderlying(AxisDirection::AxisY)] = 20;
  transformD.pretranslate(translationD);

  auto volumeD = Experimental::DetectorVolumeFactory::construct(
      portalGenerator, tContext, "VolumeD", transformD, box,
      Experimental::tryAllPortals());

  volumes = {volumeC, volumeD};
  auto containerCD = Experimental::detail::CuboidalDetectorHelper::connect(
      tContext, volumes, AxisDirection::AxisX, {}, Logging::VERBOSE);

  auto containerABCD = Experimental::detail::CuboidalDetectorHelper::connect(
      tContext, {containerAB, containerCD}, AxisDirection::AxisY, {},
      Logging::VERBOSE);

  // Check the container is constructed
  BOOST_CHECK(!containerABCD.empty());

  ObjVisualization3D obj;
  GeometryView3D::drawDetectorVolume(obj, *volumeA, tContext);
  GeometryView3D::drawDetectorVolume(obj, *volumeB, tContext);
  GeometryView3D::drawDetectorVolume(obj, *volumeC, tContext);
  GeometryView3D::drawDetectorVolume(obj, *volumeD, tContext);

  obj.write("ConnectContainers_binX.obj");
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace ActsTests

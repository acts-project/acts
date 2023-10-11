// This file is part of the Acts project.
//
// Copyright (C) 2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Geometry/GeometryContext.hpp"

namespace Acts {
namespace Experimental {
class IDetectorComponentBuilder {
 public:
  IDetectorComponentBuilder() = default;
};
}  // namespace Experimental
}  // namespace Acts

using namespace Acts;
#include "Acts/Detector/ComponentBuilderProxy.hpp"
#include "Acts/Detector/detail/CylindricalComponentProxyHelper.hpp"
#include "Acts/Tests/CommonHelpers/FloatComparisons.hpp"
#include "Acts/Utilities/Logger.hpp"

GeometryContext tContext;

BOOST_AUTO_TEST_SUITE(Detector)

BOOST_AUTO_TEST_CASE(CylindricalComponentProxyHelperTest_Z_Untoched) {
  auto idcb = std::make_shared<Experimental::IDetectorComponentBuilder>();

  auto rootProxy = Experimental::ComponentBuilderProxy::createRootProxy(
      "RootBuilder", Transform3::Identity(),
      VolumeBounds::BoundsType::eCylinder, {0.0, 10., 100.}, {binZ});

  auto child0 = rootProxy->addChildVolumeProxy(
      "Child0", Transform3::Identity() * Translation3(0., 0., -50.),
      VolumeBounds::BoundsType::eCylinder, {0.0, 10., 50.}, idcb);

  auto child1 = rootProxy->addChildVolumeProxy(
      "Child1", Transform3::Identity() * Translation3(0., 0., 50.),
      VolumeBounds::BoundsType::eCylinder, {0.0, 10., 50.}, idcb);

  Experimental::detail::CylindricalComponentProxyHelper::addGapProxies(
      *rootProxy, Acts::Logging::VERBOSE);

  Experimental::ComponentBuilderProxy::ContainerProxy container =
      std::get<Experimental::ComponentBuilderProxy::ContainerProxy>(
          rootProxy->holder);

  BOOST_CHECK(container.children.size() == 2u);
}

BOOST_AUTO_TEST_CASE(CylindricalComponentProxyHelperTest_Z_OneSide) {
  auto idcb = std::make_shared<Experimental::IDetectorComponentBuilder>();

  auto rootProxy = Experimental::ComponentBuilderProxy::createRootProxy(
      "RootBuilder", Transform3::Identity(),
      VolumeBounds::BoundsType::eCylinder, {0.0, 10., 100.}, {binZ});

  auto child0 = rootProxy->addChildVolumeProxy(
      "Child0", Transform3::Identity() * Translation3(0., 0., -50.),
      VolumeBounds::BoundsType::eCylinder, {0.0, 10., 50.}, idcb);

  auto child1 = rootProxy->addChildVolumeProxy(
      "Child1", Transform3::Identity() * Translation3(0., 0., 50.),
      VolumeBounds::BoundsType::eCylinder, {0.0, 10., 30.}, idcb);

  Experimental::detail::CylindricalComponentProxyHelper::addGapProxies(
      *rootProxy, Acts::Logging::VERBOSE);

  Experimental::ComponentBuilderProxy::ContainerProxy container =
      std::get<Experimental::ComponentBuilderProxy::ContainerProxy>(
          rootProxy->holder);

  BOOST_CHECK(container.children.size() == 4u);
}

BOOST_AUTO_TEST_CASE(CylindricalComponentProxyHelperTest_Z_BothSides) {
  auto idcb = std::make_shared<Experimental::IDetectorComponentBuilder>();

  auto rootProxy = Experimental::ComponentBuilderProxy::createRootProxy(
      "RootBuilder", Transform3::Identity(),
      VolumeBounds::BoundsType::eCylinder, {0.0, 10., 100.}, {binZ});

  auto child0 = rootProxy->addChildVolumeProxy(
      "Child0", Transform3::Identity() * Translation3(0., 0., -50.),
      VolumeBounds::BoundsType::eCylinder, {0.0, 10., 20.}, idcb);

  auto child1 = rootProxy->addChildVolumeProxy(
      "Child1", Transform3::Identity() * Translation3(0., 0., 50.),
      VolumeBounds::BoundsType::eCylinder, {0.0, 10., 20.}, idcb);

  Experimental::detail::CylindricalComponentProxyHelper::addGapProxies(
      *rootProxy, Acts::Logging::VERBOSE);

  Experimental::ComponentBuilderProxy::ContainerProxy container =
      std::get<Experimental::ComponentBuilderProxy::ContainerProxy>(
          rootProxy->holder);

  BOOST_CHECK(container.children.size() == 5u);

  auto gapProxy0 = container.children[0];
  auto childProxy0 = container.children[1];
  auto gapProxy1 = container.children[2];
  auto childProxy1 = container.children[3];
  auto gapProxy2 = container.children[4];

  // Gaps are positioned
  CHECK_CLOSE_ABS(gapProxy0->transform.translation().z(), -85., 1e-6);
  CHECK_CLOSE_ABS(childProxy0->transform.translation().z(), -50., 1e-6);
  CHECK_CLOSE_ABS(gapProxy1->transform.translation().z(), 0., 1e-6);
  CHECK_CLOSE_ABS(childProxy1->transform.translation().z(), 50., 1e-6);
  CHECK_CLOSE_ABS(gapProxy2->transform.translation().z(), 85., 1e-6);

  // Gaps have the right length
  CHECK_CLOSE_ABS(gapProxy0->boundaryValues[2], 15., 1e-6);
  CHECK_CLOSE_ABS(gapProxy1->boundaryValues[2], 30., 1e-6);
  CHECK_CLOSE_ABS(gapProxy2->boundaryValues[2], 15., 1e-6);
}

BOOST_AUTO_TEST_CASE(CylindricalComponentProxyHelperTest_R_Untouched) {
  auto idcb = std::make_shared<Experimental::IDetectorComponentBuilder>();

  auto rootProxy = Experimental::ComponentBuilderProxy::createRootProxy(
      "RootBuilder", Transform3::Identity(),
      VolumeBounds::BoundsType::eCylinder, {0.0, 50., 100.}, {binR});

  auto child0 = rootProxy->addChildVolumeProxy(
      "Child0", Transform3::Identity(), VolumeBounds::BoundsType::eCylinder,
      {0.0, 10., 100.}, idcb);

  auto child1 = rootProxy->addChildVolumeProxy(
      "Child1", Transform3::Identity(), VolumeBounds::BoundsType::eCylinder,
      {10., 50., 100.}, idcb);

  Experimental::detail::CylindricalComponentProxyHelper::addGapProxies(
      *rootProxy, Acts::Logging::VERBOSE);

  Experimental::ComponentBuilderProxy::ContainerProxy container =
      std::get<Experimental::ComponentBuilderProxy::ContainerProxy>(
          rootProxy->holder);

  BOOST_CHECK(container.children.size() == 2u);
}

BOOST_AUTO_TEST_CASE(CylindricalComponentProxyHelperTest_R_Inner) {
  auto idcb = std::make_shared<Experimental::IDetectorComponentBuilder>();

  auto rootProxy = Experimental::ComponentBuilderProxy::createRootProxy(
      "RootBuilder", Transform3::Identity(),
      VolumeBounds::BoundsType::eCylinder, {0.0, 50., 100.}, {binR});

  auto child0 = rootProxy->addChildVolumeProxy(
      "Child0", Transform3::Identity(), VolumeBounds::BoundsType::eCylinder,
      {10., 20., 100.}, idcb);

  auto child1 = rootProxy->addChildVolumeProxy(
      "Child1", Transform3::Identity(), VolumeBounds::BoundsType::eCylinder,
      {20., 50., 100.}, idcb);

  Experimental::detail::CylindricalComponentProxyHelper::addGapProxies(
      *rootProxy, Acts::Logging::VERBOSE);

  Experimental::ComponentBuilderProxy::ContainerProxy container =
      std::get<Experimental::ComponentBuilderProxy::ContainerProxy>(
          rootProxy->holder);

  BOOST_CHECK(container.children.size() == 3u);
  auto gapProxy0 = container.children[0];

  CHECK_CLOSE_ABS(gapProxy0->boundaryValues[0], 0., 1e-6);
  CHECK_CLOSE_ABS(gapProxy0->boundaryValues[1], 10., 1e-6);
}

BOOST_AUTO_TEST_CASE(CylindricalComponentProxyHelperTest_R_All) {
  auto idcb = std::make_shared<Experimental::IDetectorComponentBuilder>();

  auto rootProxy = Experimental::ComponentBuilderProxy::createRootProxy(
      "RootBuilder", Transform3::Identity(),
      VolumeBounds::BoundsType::eCylinder, {0.0, 50., 100.}, {binR});

  auto child0 = rootProxy->addChildVolumeProxy(
      "Child0", Transform3::Identity(), VolumeBounds::BoundsType::eCylinder,
      {10., 20., 100.}, idcb);

  auto child1 = rootProxy->addChildVolumeProxy(
      "Child1", Transform3::Identity(), VolumeBounds::BoundsType::eCylinder,
      {40., 45., 100.}, idcb);

  Experimental::detail::CylindricalComponentProxyHelper::addGapProxies(
      *rootProxy, Acts::Logging::VERBOSE);

  Experimental::ComponentBuilderProxy::ContainerProxy container =
      std::get<Experimental::ComponentBuilderProxy::ContainerProxy>(
          rootProxy->holder);

  BOOST_CHECK(container.children.size() == 5u);
  auto gapProxy0 = container.children[0];
  auto gapProxy1 = container.children[2];
  auto gapProxy2 = container.children[4];

  CHECK_CLOSE_ABS(gapProxy0->boundaryValues[0], 0., 1e-6);
  CHECK_CLOSE_ABS(gapProxy0->boundaryValues[1], 10., 1e-6);
  CHECK_CLOSE_ABS(gapProxy1->boundaryValues[0], 20., 1e-6);
  CHECK_CLOSE_ABS(gapProxy1->boundaryValues[1], 40., 1e-6);
  CHECK_CLOSE_ABS(gapProxy2->boundaryValues[0], 45., 1e-6);
  CHECK_CLOSE_ABS(gapProxy2->boundaryValues[1], 50., 1e-6);
}

BOOST_AUTO_TEST_SUITE_END()

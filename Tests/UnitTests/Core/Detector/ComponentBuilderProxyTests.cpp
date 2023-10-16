// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Detector/ComponentBuilderProxy.hpp"
#include "Acts/Geometry/VolumeBounds.hpp"

#include <exception>
#include <memory>
#include <variant>

namespace Acts {
namespace Experimental {
class IDetectorComponentBuilder {
 public:
  IDetectorComponentBuilder() = default;
};
}  // namespace Experimental
}  // namespace Acts

using namespace Acts;
using namespace Acts::Experimental;

BOOST_AUTO_TEST_SUITE(Experimental)

BOOST_AUTO_TEST_CASE(ComponentBuilderProxyTest) {
  auto idcb = std::make_shared<IDetectorComponentBuilder>();

  auto singleProxy = ComponentBuilderProxy::createRootProxy(
      "SingleProxy", Transform3::Identity(),
      VolumeBounds::BoundsType::eCylinder, {0.0, 10., 100.}, {}, idcb);

  auto rootProxy = ComponentBuilderProxy::createRootProxy(
      "RootBuilder", Transform3::Identity(),
      VolumeBounds::BoundsType::eCylinder, {0.0, 10., 100.}, {binZ});

  // The root proxy has no parent
  BOOST_CHECK(rootProxy->parent == nullptr);

  // Works - add a child to the root proxy
  rootProxy->addChildVolumeProxy(
      "FirstChild", Transform3::Identity() * Translation3(0., 0., -50.),
      VolumeBounds::BoundsType::eCylinder, {0.0, 10., 50.}, idcb);

  // Does not work - add a child to the single proxy
  BOOST_CHECK_THROW(
      singleProxy->addChildVolumeProxy(
          "FirstChild", Transform3::Identity() * Translation3(0., 0., -50.),
          VolumeBounds::BoundsType::eCylinder, {0.0, 10., 50.}, idcb),
      std::runtime_error);

  // The root proxy now turned into an container proxy
  BOOST_CHECK(std::holds_alternative<ComponentBuilderProxy::ContainerProxy>(
      rootProxy->holder));

  // Works - add a child container to the root proxy
  rootProxy->addChildContainerProxy(
      "SecondChild", Transform3::Identity() * Translation3(0., 0., 50),
      VolumeBounds::BoundsType::eCylinder, {0.0, 10., 50.});

  auto& cProxy =
      std::get<ComponentBuilderProxy::ContainerProxy>(rootProxy->holder);
  BOOST_CHECK(cProxy.children.size() == 2);

  // Does not work - add a child container to the single proxy
  BOOST_CHECK_THROW(
      singleProxy->addChildContainerProxy(
          "SecondChild", Transform3::Identity() * Translation3(0., 0., 50),
          VolumeBounds::BoundsType::eCylinder, {0.0, 10., 50.}),
      std::runtime_error);
}

BOOST_AUTO_TEST_SUITE_END()

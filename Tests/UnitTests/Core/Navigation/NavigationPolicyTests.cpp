// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Definitions/Units.hpp"
#include "Acts/Geometry/CylinderVolumeBounds.hpp"
#include "Acts/Geometry/NavigationPolicyFactory.hpp"
#include "Acts/Geometry/TrackingVolume.hpp"
#include "Acts/Navigation/INavigationPolicy.hpp"
#include "Acts/Navigation/MultiNavigationPolicy.hpp"
#include "Acts/Navigation/NavigationDelegate.hpp"
#include "Acts/Navigation/NavigationStream.hpp"
#include "Acts/Navigation/TryAllNavigationPolicies.hpp"

using namespace Acts;
using namespace Acts::UnitLiterals;

BOOST_AUTO_TEST_SUITE(NavigationPolicyTests)

struct APolicy : public INavigationPolicy {
  APolicy(const TrackingVolume& /*volume*/) {}

  void updateState(const NavigationArguments& /*unused*/) const {
    const_cast<APolicy*>(this)->executed = true;
  }

  void connect(NavigationDelegate& delegate) const override {
    connectDefault<TryAllPortalNavigationPolicy>(delegate);
  }

  bool executed = false;
};

struct BPolicy : public INavigationPolicy {
  struct Config {
    int value;
  };

  BPolicy(const TrackingVolume& /*volume*/, Config config) : m_config(config) {}

  void connect(NavigationDelegate& delegate) const override {
    connectDefault<TryAllPortalNavigationPolicy>(delegate);
  }

  void updateState(const NavigationArguments& /*unused*/) const {
    const_cast<BPolicy*>(this)->executed = true;
    const_cast<BPolicy*>(this)->value = m_config.value;
  }

  bool executed = false;
  int value = 0;

  Config m_config;
};

BOOST_AUTO_TEST_CASE(DirectTest) {
  TrackingVolume volume{
      Transform3::Identity(),
      std::make_shared<CylinderVolumeBounds>(250_mm, 400_mm, 310_mm),
      "PixelLayer3"};

  MultiNavigationPolicy policy{APolicy{volume},
                               BPolicy{volume, {.value = 4242}}};

  NavigationDelegate delegate;
  policy.connect(delegate);

  NavigationStream main;
  delegate(NavigationArguments{
      .main = main, .position = Vector3::Zero(), .direction = Vector3::Zero()});

  BOOST_CHECK(std::get<APolicy>(policy.policies()).executed);
  BOOST_CHECK(std::get<BPolicy>(policy.policies()).executed);
  BOOST_CHECK_EQUAL(std::get<BPolicy>(policy.policies()).value, 4242);
}

BOOST_AUTO_TEST_CASE(FactoryTest) {
  TrackingVolume volume{
      Transform3::Identity(),
      std::make_shared<CylinderVolumeBounds>(250_mm, 400_mm, 310_mm),
      "PixelLayer3"};

  BPolicy::Config config{.value = 42};

  std::function<std::unique_ptr<INavigationPolicy>(const TrackingVolume&)>
      factory = NavigationPolicyFactory::make()
                    .add<APolicy>()         // no arguments
                    .add<BPolicy>(config);  // config struct as argument

  auto policyBase = factory(volume);
  auto policyBase2 = factory(volume);

  auto& policy =
      dynamic_cast<MultiNavigationPolicy<APolicy, BPolicy>&>(*policyBase);

  NavigationDelegate delegate;
  policy.connect(delegate);

  NavigationStream main;
  delegate(NavigationArguments{
      .main = main, .position = Vector3::Zero(), .direction = Vector3::Zero()});

  BOOST_CHECK(std::get<APolicy>(policy.policies()).executed);
  BOOST_CHECK(std::get<BPolicy>(policy.policies()).executed);
  BOOST_CHECK_EQUAL(std::get<BPolicy>(policy.policies()).value, 42);

  auto& policy2 =
      dynamic_cast<MultiNavigationPolicy<APolicy, BPolicy>&>(*policyBase2);

  NavigationDelegate delegate2;
  policyBase2->connect(delegate2);

  delegate2(NavigationArguments{
      .main = main, .position = Vector3::Zero(), .direction = Vector3::Zero()});

  BOOST_CHECK(std::get<APolicy>(policy2.policies()).executed);
  BOOST_CHECK(std::get<BPolicy>(policy2.policies()).executed);
  BOOST_CHECK_EQUAL(std::get<BPolicy>(policy2.policies()).value, 42);
}

BOOST_AUTO_TEST_CASE(AsUniquePtrTest) {
  TrackingVolume volume{
      Transform3::Identity(),
      std::make_shared<CylinderVolumeBounds>(250_mm, 400_mm, 310_mm),
      "PixelLayer3"};

  std::unique_ptr<NavigationPolicyFactory> factory =
      NavigationPolicyFactory::make().add<APolicy>().asUniquePtr();

  auto policyBase = factory->build(volume);
  auto& policy = dynamic_cast<MultiNavigationPolicy<APolicy>&>(*policyBase);

  NavigationDelegate delegate;
  policyBase->connect(delegate);

  NavigationStream main;
  delegate(NavigationArguments{
      .main = main, .position = Vector3::Zero(), .direction = Vector3::Zero()});

  BOOST_CHECK(std::get<APolicy>(policy.policies()).executed);
}

BOOST_AUTO_TEST_SUITE_END()

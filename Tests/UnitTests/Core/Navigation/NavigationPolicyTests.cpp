// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/Tolerance.hpp"
#include "Acts/Definitions/Units.hpp"
#include "Acts/Geometry/CylinderPortalShell.hpp"
#include "Acts/Geometry/CylinderVolumeBounds.hpp"
#include "Acts/Geometry/NavigationPolicyFactory.hpp"
#include "Acts/Geometry/TrackingVolume.hpp"
#include "Acts/Navigation/CylinderNavigationPolicy.hpp"
#include "Acts/Navigation/INavigationPolicy.hpp"
#include "Acts/Navigation/MultiNavigationPolicy.hpp"
#include "Acts/Navigation/NavigationDelegate.hpp"
#include "Acts/Navigation/NavigationStream.hpp"
#include "Acts/Navigation/TryAllNavigationPolicy.hpp"

using namespace Acts;
using namespace Acts::UnitLiterals;

BOOST_AUTO_TEST_SUITE(NavigationPolicyTests)

GeometryContext gctx;
auto logger = getDefaultLogger("NavigationPolicyTests", Logging::VERBOSE);

struct APolicy : public INavigationPolicy {
  APolicy(const GeometryContext& /*gctx*/, const TrackingVolume& /*volume*/,
          const Logger& /*logger*/) {}

  void initializeCandidates(const NavigationArguments& /*unused*/,
                            AppendOnlyNavigationStream& /*unused*/,
                            const Logger& /*unused*/) const {
    const_cast<APolicy*>(this)->executed = true;
  }

  void connect(NavigationDelegate& delegate) const override {
    connectDefault<APolicy>(delegate);
  }

  bool executed = false;
};

struct BPolicy : public INavigationPolicy {
  struct Config {
    int value;
  };

  BPolicy(const GeometryContext& /*gctx*/, const TrackingVolume& /*volume*/,
          const Logger& /*logger*/, Config config)
      : m_config(config) {}

  void connect(NavigationDelegate& delegate) const override {
    connectDefault<BPolicy>(delegate);
  }

  void initializeCandidates(const NavigationArguments& /*unused*/,
                            AppendOnlyNavigationStream& /*unused*/,
                            const Logger& /*unused*/) const {
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

  MultiNavigationPolicy policy{APolicy{gctx, volume, *logger},
                               BPolicy{gctx, volume, *logger, {.value = 4242}}};

  NavigationDelegate delegate;
  policy.connect(delegate);

  NavigationStream main;
  AppendOnlyNavigationStream stream{main};
  delegate(NavigationArguments{.position = Vector3::Zero(),
                               .direction = Vector3::Zero()},
           stream, *logger);

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

  std::function<std::unique_ptr<INavigationPolicy>(
      const GeometryContext&, const TrackingVolume&, const Logger&)>
      factory = NavigationPolicyFactory::make()
                    .add<APolicy>()         // no arguments
                    .add<BPolicy>(config);  // config struct as argument

  auto policyBase = factory(gctx, volume, *logger);
  auto policyBase2 = factory(gctx, volume, *logger);

  auto& policy =
      dynamic_cast<MultiNavigationPolicy<APolicy, BPolicy>&>(*policyBase);

  NavigationDelegate delegate;
  policy.connect(delegate);

  NavigationStream main;
  AppendOnlyNavigationStream stream{main};
  delegate(NavigationArguments{.position = Vector3::Zero(),
                               .direction = Vector3::Zero()},
           stream, *logger);

  BOOST_CHECK(std::get<APolicy>(policy.policies()).executed);
  BOOST_CHECK(std::get<BPolicy>(policy.policies()).executed);
  BOOST_CHECK_EQUAL(std::get<BPolicy>(policy.policies()).value, 42);

  auto& policy2 =
      dynamic_cast<MultiNavigationPolicy<APolicy, BPolicy>&>(*policyBase2);

  NavigationDelegate delegate2;
  policyBase2->connect(delegate2);

  delegate2(NavigationArguments{.position = Vector3::Zero(),
                                .direction = Vector3::Zero()},
            stream, *logger);

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

  auto policyBase = factory->build(gctx, volume, *logger);
  auto& policy = dynamic_cast<MultiNavigationPolicy<APolicy>&>(*policyBase);

  NavigationDelegate delegate;
  policyBase->connect(delegate);

  NavigationStream main;
  AppendOnlyNavigationStream stream{main};
  delegate(NavigationArguments{.position = Vector3::Zero(),
                               .direction = Vector3::Zero()},
           stream, *logger);

  BOOST_CHECK(std::get<APolicy>(policy.policies()).executed);
}

struct CPolicy : public INavigationPolicy {};

template <typename T>
struct CPolicySpecialized : public CPolicy {
  struct Config {
    T value;
  };

  CPolicySpecialized(const TrackingVolume& /*volume*/, Config config)
      : m_config(config) {}

  void connect(NavigationDelegate& delegate) const override {
    connectDefault<CPolicySpecialized<T>>(delegate);
  }

  void initializeCandidates(const NavigationArguments& /*unused*/,
                            AppendOnlyNavigationStream& /*stream*/,
                            const Logger& /*logger*/) const {
    auto* self = const_cast<CPolicySpecialized<int>*>(this);
    self->executed = true;
    self->value = m_config.value;
  }

  bool executed = false;
  int value = 0;

  Config m_config;
};

struct IsolatedConfig {
  int value;
};

auto makeCPolicy(const GeometryContext& /*gctx*/, const TrackingVolume& volume,
                 const Logger& /*logger*/, IsolatedConfig config) {
  // I can do arbitrary stuff here
  CPolicySpecialized<int>::Config config2{.value = config.value};
  return CPolicySpecialized<int>(volume, config2);
}

BOOST_AUTO_TEST_CASE(IsolatedFactory) {
  TrackingVolume volume{
      Transform3::Identity(),
      std::make_shared<CylinderVolumeBounds>(250_mm, 400_mm, 310_mm),
      "PixelLayer3"};

  IsolatedConfig config{.value = 44};
  auto factory =
      NavigationPolicyFactory::make().add<APolicy>().add(makeCPolicy, config);

  auto factory2 =
      NavigationPolicyFactory::make().add(makeCPolicy, config).add<APolicy>();

  auto policyBase = factory(gctx, volume, *logger);
  auto& policy =
      dynamic_cast<MultiNavigationPolicy<APolicy, CPolicySpecialized<int>>&>(
          *policyBase);

  NavigationDelegate delegate;
  policyBase->connect(delegate);

  NavigationStream main;
  AppendOnlyNavigationStream stream{main};
  delegate(NavigationArguments{.position = Vector3::Zero(),
                               .direction = Vector3::Zero()},
           stream, *logger);

  BOOST_CHECK(std::get<APolicy>(policy.policies()).executed);
  BOOST_CHECK(std::get<CPolicySpecialized<int>>(policy.policies()).executed);
  BOOST_CHECK_EQUAL(std::get<CPolicySpecialized<int>>(policy.policies()).value,
                    44);
}

BOOST_AUTO_TEST_CASE(CylinderPolicyTest) {
  using enum CylinderVolumeBounds::Face;

  Transform3 transform = Transform3::Identity();
  // transform *= AngleAxis3{std::numbers::pi / 2, Vector3::UnitY()};
  auto cylBounds =
      std::make_shared<CylinderVolumeBounds>(100_mm, 400_mm, 300_mm);
  auto cylVolume =
      std::make_shared<TrackingVolume>(transform, cylBounds, "CylinderVolume");
  SingleCylinderPortalShell shell{*cylVolume};
  shell.applyToVolume();

  auto getTruth = [&](const Vector3& position, const Vector3& direction,
                      bool posOnly = true) {
    TryAllNavigationPolicy tryAll(gctx, *cylVolume, *logger);
    NavigationArguments args{.position = position, .direction = direction};
    NavigationStream main;
    AppendOnlyNavigationStream stream{main};
    tryAll.initializeCandidates(args, stream, *logger);
    main.initialize(gctx, {position, direction}, BoundaryTolerance::None());
    std::vector<const Portal*> portals;
    for (auto& candidate : main.candidates()) {
      if (!candidate.intersection.isValid()) {
        continue;
      }

      if (posOnly &&
          !detail::checkPathLength(
              candidate.intersection.pathLength(), s_onSurfaceTolerance,
              std::numeric_limits<double>::max(), *logger)) {
        continue;
      }

      portals.push_back(candidate.portal);
    }

    auto outerIt = std::ranges::find(portals, shell.portal(OuterCylinder));
    auto innerIt = std::ranges::find(portals, shell.portal(InnerCylinder));

    if (outerIt != portals.end() && innerIt != portals.end()) {
      // have both inner and outer cylinder. we only expect the inner one to be
      // returned
      portals.erase(outerIt);
    }

    return portals;
  };

  CylinderNavigationPolicy policy(gctx, *cylVolume, *logger, {});
  auto getSmart = [&](const Vector3& position, const Vector3& direction) {
    NavigationArguments args{.position = position, .direction = direction};
    NavigationStream main;
    AppendOnlyNavigationStream stream{main};
    policy.initializeCandidates(args, stream, *logger);

    std::vector<const Portal*> portals;
    // We don't filter here, because we want to test the candidates as they come
    // out of the policy
    for (auto& candidate : main.candidates()) {
      portals.push_back(candidate.portal);
    }
    return portals;
  };

  auto checkEqual = [&](const std::vector<const Portal*>& exp,
                        const std::vector<const Portal*>& act) {
    auto which = [&](const Portal* p) -> std::string {
      if (p == shell.portal(InnerCylinder)) {
        return "InnerCylinder";
      }
      if (p == shell.portal(OuterCylinder)) {
        return "OuterCylinder";
      }
      if (p == shell.portal(PositiveDisc)) {
        return "PositiveDisc";
      }
      if (p == shell.portal(NegativeDisc)) {
        return "NegativeDisc";
      }
      BOOST_FAIL("Unknown portal");
      return "";  // unreachable
    };

    std::size_t imin = std::min(exp.size(), act.size());
    for (std::size_t i = 0; i < imin; ++i) {
      const auto& p1 = exp.at(i);
      const auto& p2 = act.at(i);
      if (p1 != p2) {
        BOOST_ERROR("Portal mismatch at index " << i << " exp: " << which(p1)
                                                << " act: " << which(p2));
      }
    }

    if (exp.size() > imin) {
      BOOST_ERROR("Expected more portals than actual");
      for (std::size_t i = imin; i < exp.size(); ++i) {
        BOOST_ERROR("Expected additional portal " << which(exp.at(i)));
      }
    }

    if (act.size() > imin) {
      BOOST_ERROR("Actual more portals than expected");
      for (std::size_t i = imin; i < act.size(); ++i) {
        BOOST_ERROR("Actual additional portal " << which(act.at(i)));
      }
    }
  };

  {
    Vector3 position = Vector3::UnitX() * 150_mm;
    Vector3 direction = Vector3::UnitZ();

    auto exp = getTruth(position, direction);

    BOOST_CHECK(exp.size() == 1);
    BOOST_CHECK(exp.at(0) == shell.portal(PositiveDisc));

    auto act = getSmart(position, direction);
    checkEqual(exp, act);
  }

  {
    Vector3 position = Vector3::UnitX() * 150_mm;
    Vector3 direction = Vector3{1, 1, 0}.normalized();

    auto exp = getTruth(position, direction);

    BOOST_CHECK(exp.size() == 1);
    BOOST_CHECK(exp.at(0) == shell.portal(OuterCylinder));

    auto act = getSmart(position, direction);
    checkEqual(exp, act);
  }

  {
    Vector3 position = Vector3::UnitX() * 150_mm;
    Vector3 direction = Vector3{-1, 0, 0}.normalized();

    auto exp = getTruth(position, direction);

    BOOST_CHECK(exp.size() == 1);
    BOOST_CHECK(exp.at(0) == shell.portal(InnerCylinder));

    auto act = getSmart(position, direction);
    checkEqual(exp, act);
  }

  {
    Vector3 position = Vector3::UnitX() * 150_mm;
    Vector3 direction = -Vector3::UnitZ();

    auto exp = getTruth(position, direction);

    BOOST_CHECK(exp.size() == 1);
    BOOST_CHECK(exp.at(0) == shell.portal(NegativeDisc));

    auto act = getSmart(position, direction);
    checkEqual(exp, act);
  }

  {
    Vector3 position{50, -200, 0};
    Vector3 direction = Vector3{0, 1.5, 1}.normalized();

    auto exp = getTruth(position, direction);

    BOOST_CHECK(exp.size() == 2);
    BOOST_CHECK(exp.at(0) == shell.portal(InnerCylinder));
    BOOST_CHECK(exp.at(1) == shell.portal(PositiveDisc));

    auto act = getSmart(position, direction);
    checkEqual(exp, act);
  }

  {
    Vector3 position{50, -200, 0};
    Vector3 direction = Vector3{0, 1.2, 1}.normalized();

    auto exp = getTruth(position, direction);

    BOOST_CHECK(exp.size() == 2);
    BOOST_CHECK(exp.at(0) == shell.portal(InnerCylinder));
    BOOST_CHECK(exp.at(1) == shell.portal(PositiveDisc));

    auto act = getSmart(position, direction);
    checkEqual(exp, act);
  }

  {
    Vector3 position{50, -200, 0};
    Vector3 direction = Vector3{0, 0.9, 1}.normalized();

    auto exp = getTruth(position, direction);

    BOOST_CHECK(exp.size() == 1);
    BOOST_CHECK(exp.at(0) == shell.portal(InnerCylinder));

    auto act = getSmart(position, direction);
    checkEqual(exp, act);
  }

  {
    Vector3 position{20, -200, 0};
    Vector3 direction = Vector3{0.45, 0.9, 1}.normalized();

    auto exp = getTruth(position, direction);

    BOOST_CHECK(exp.size() == 1);
    BOOST_CHECK(exp.at(0) == shell.portal(PositiveDisc));

    auto act = getSmart(position, direction);
    checkEqual(exp, act);
  }

  {
    Vector3 position{20, -200, 0};
    Vector3 direction = Vector3{0.45, 0.9, -1}.normalized();

    auto exp = getTruth(position, direction);

    BOOST_CHECK(exp.size() == 1);
    BOOST_CHECK(exp.at(0) == shell.portal(NegativeDisc));

    auto act = getSmart(position, direction);
    checkEqual(exp, act);
  }

  {
    Vector3 position{400 * std::cos(std::numbers::pi / 4),
                     400 * std::sin(std::numbers::pi / 4), 0};
    Vector3 direction = Vector3{0.45, -0.9, -0.1}.normalized();

    // We're sitting ON the outer cylinder here and missing the disc and inner
    // cylinder: need to return the outer cylinder
    auto exp = getTruth(position, direction, false);

    BOOST_CHECK(exp.size() == 1);
    BOOST_CHECK(exp.at(0) == shell.portal(OuterCylinder));

    auto act = getSmart(position, direction);
    checkEqual(exp, act);
  }

  {
    Vector3 position{400 * std::cos(std::numbers::pi / 4),
                     400 * std::sin(std::numbers::pi / 4), 0};
    Vector3 direction = Vector3{-0.3, -0.9, -0.1}.normalized();

    // We're sitting ON the outer cylinder here and missing the disc and inner
    // cylinder: need to return the outer cylinder
    auto exp = getTruth(position, direction, false);

    BOOST_CHECK(exp.size() == 1);
    BOOST_CHECK(exp.at(0) == shell.portal(OuterCylinder));

    auto act = getSmart(position, direction);
    checkEqual(exp, act);
  }

  {
    Vector3 position{100 * std::cos(std::numbers::pi / 4),
                     100 * std::sin(std::numbers::pi / 4), 0};
    double dangle = 0.1;
    Vector3 direction = Vector3{std::cos(dangle), std::sin(dangle), 0.01};

    // We're sitting ON the Inner cylinder here and  pointing outwards
    auto exp = getTruth(position, direction);

    BOOST_CHECK(exp.size() == 1);
    BOOST_CHECK(exp.at(0) == shell.portal(OuterCylinder));

    auto act = getSmart(position, direction);
    checkEqual(exp, act);
  }
}

BOOST_AUTO_TEST_SUITE_END()

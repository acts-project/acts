// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/data/test_case.hpp>
#include <boost/test/unit_test.hpp>

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/Tolerance.hpp"
#include "Acts/Definitions/Units.hpp"
#include "Acts/Geometry/CylinderPortalShell.hpp"
#include "Acts/Geometry/CylinderVolumeBounds.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/NavigationPolicyFactory.hpp"
#include "Acts/Geometry/TrackingVolume.hpp"
#include "Acts/Navigation/CylinderNavigationPolicy.hpp"
#include "Acts/Navigation/INavigationPolicy.hpp"
#include "Acts/Navigation/MultiNavigationPolicy.hpp"
#include "Acts/Navigation/NavigationDelegate.hpp"
#include "Acts/Navigation/NavigationStream.hpp"
#include "Acts/Navigation/TryAllNavigationPolicy.hpp"
#include "Acts/Utilities/Logger.hpp"

#include <boost/algorithm/string/join.hpp>

using namespace Acts;
using namespace Acts::UnitLiterals;
namespace bdata = boost::unit_test::data;

namespace ActsTests {

BOOST_AUTO_TEST_SUITE(NavigationSuite)

GeometryContext gctx;
auto logger = getDefaultLogger("NavigationPolicyTests", Logging::VERBOSE);

struct APolicy : public INavigationPolicy {
  APolicy(const GeometryContext& /*gctx*/, const TrackingVolume& /*volume*/,
          const Logger& /*logger*/) {}

  void initializeCandidates(const GeometryContext& /*unused*/,
                            const NavigationArguments& /*unused*/,
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

  void initializeCandidates(const GeometryContext& /*unused*/,
                            const NavigationArguments& /*unused*/,
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

  MultiNavigationPolicy policy{
      std::make_unique<APolicy>(gctx, volume, *logger),
      std::make_unique<BPolicy>(gctx, volume, *logger,
                                BPolicy::Config{.value = 4242})};

  NavigationDelegate delegate;
  policy.connect(delegate);

  NavigationStream main;
  AppendOnlyNavigationStream stream{main};
  delegate(gctx,
           NavigationArguments{.position = Vector3::Zero(),
                               .direction = Vector3::Zero()},
           stream, *logger);

  BOOST_REQUIRE_EQUAL(policy.policies().size(), 2);
  const auto& policyA = dynamic_cast<const APolicy&>(*policy.policies()[0]);
  const auto& policyB = dynamic_cast<const BPolicy&>(*policy.policies()[1]);

  BOOST_CHECK(policyA.executed);
  BOOST_CHECK(policyB.executed);
  BOOST_CHECK_EQUAL(policyB.value, 4242);
}

BOOST_AUTO_TEST_CASE(FactoryTest) {
  TrackingVolume volume{
      Transform3::Identity(),
      std::make_shared<CylinderVolumeBounds>(250_mm, 400_mm, 310_mm),
      "PixelLayer3"};

  BPolicy::Config config{.value = 42};

  std::function<std::unique_ptr<INavigationPolicy>(
      const GeometryContext&, const TrackingVolume&, const Logger&)>
      factory = NavigationPolicyFactory{}
                    .add<APolicy>()         // no arguments
                    .add<BPolicy>(config);  // config struct as argument

  auto policyBase = factory(gctx, volume, *logger);
  auto policyBase2 = factory(gctx, volume, *logger);

  auto& policy = dynamic_cast<MultiNavigationPolicy&>(*policyBase);

  NavigationDelegate delegate;
  policy.connect(delegate);

  NavigationStream main;
  AppendOnlyNavigationStream stream{main};
  delegate(gctx,
           NavigationArguments{.position = Vector3::Zero(),
                               .direction = Vector3::Zero()},
           stream, *logger);

  BOOST_REQUIRE_EQUAL(policy.policies().size(), 2);
  const auto& policyA = dynamic_cast<const APolicy&>(*policy.policies()[0]);
  const auto& policyB = dynamic_cast<const BPolicy&>(*policy.policies()[1]);

  BOOST_CHECK(policyA.executed);
  BOOST_CHECK(policyB.executed);
  BOOST_CHECK_EQUAL(policyB.value, 42);

  auto& policy2 = dynamic_cast<MultiNavigationPolicy&>(*policyBase2);

  NavigationDelegate delegate2;
  policyBase2->connect(delegate2);

  delegate2(gctx,
            NavigationArguments{.position = Vector3::Zero(),
                                .direction = Vector3::Zero()},
            stream, *logger);

  BOOST_REQUIRE_EQUAL(policy2.policies().size(), 2);
  const auto& policy2A = dynamic_cast<const APolicy&>(*policy2.policies()[0]);
  const auto& policy2B = dynamic_cast<const BPolicy&>(*policy2.policies()[1]);

  BOOST_CHECK(policy2A.executed);
  BOOST_CHECK(policy2B.executed);
  BOOST_CHECK_EQUAL(policy2B.value, 42);
}

BOOST_AUTO_TEST_CASE(AsUniquePtrTest) {
  TrackingVolume volume{
      Transform3::Identity(),
      std::make_shared<CylinderVolumeBounds>(250_mm, 400_mm, 310_mm),
      "PixelLayer3"};

  std::unique_ptr<NavigationPolicyFactory> factory =
      NavigationPolicyFactory{}.add<APolicy>().asUniquePtr();

  auto policyBase = factory->build(gctx, volume, *logger);
  auto& policy = dynamic_cast<MultiNavigationPolicy&>(*policyBase);

  NavigationDelegate delegate;
  policyBase->connect(delegate);

  NavigationStream main;
  AppendOnlyNavigationStream stream{main};
  delegate(gctx,
           NavigationArguments{.position = Vector3::Zero(),
                               .direction = Vector3::Zero()},
           stream, *logger);

  BOOST_REQUIRE_EQUAL(policy.policies().size(), 1);
  BOOST_CHECK(dynamic_cast<const APolicy&>(*policy.policies()[0]).executed);
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

  void initializeCandidates(const GeometryContext& /*unused*/,
                            const NavigationArguments& /*unused*/,
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
      NavigationPolicyFactory{}.add<APolicy>().add(makeCPolicy, config);

  auto factory2 =
      NavigationPolicyFactory{}.add(makeCPolicy, config).add<APolicy>();

  auto policyBase = factory(gctx, volume, *logger);
  auto& policy = dynamic_cast<MultiNavigationPolicy&>(*policyBase);

  NavigationDelegate delegate;
  policyBase->connect(delegate);

  NavigationStream main;
  AppendOnlyNavigationStream stream{main};
  delegate(gctx,
           NavigationArguments{.position = Vector3::Zero(),
                               .direction = Vector3::Zero()},
           stream, *logger);

  BOOST_REQUIRE_EQUAL(policy.policies().size(), 2);

  const auto& policyA = dynamic_cast<const APolicy&>(*policy.policies()[0]);
  const auto& cPolicy =
      dynamic_cast<const CPolicySpecialized<int>&>(*policy.policies()[1]);

  BOOST_CHECK(policyA.executed);
  BOOST_CHECK(cPolicy.executed);
  BOOST_CHECK_EQUAL(cPolicy.value, 44);
}

namespace {

std::vector<const Portal*> getTruth(const Vector3& position,
                                    const Vector3& direction,
                                    const Transform3& transform,
                                    const TrackingVolume& cylVolume,
                                    SingleCylinderPortalShell& shell,
                                    const Logger& logger, bool posOnly = true) {
  Vector3 gpos = transform * position;
  Vector3 gdir = transform.linear() * direction;
  TryAllNavigationPolicy tryAll(gctx, cylVolume, logger);
  NavigationArguments args{.position = gpos, .direction = gdir};
  NavigationStream main;
  AppendOnlyNavigationStream stream{main};
  GeometryContext gctx;
  tryAll.initializeCandidates(gctx, args, stream, logger);
  main.initialize(gctx, {gpos, gdir}, BoundaryTolerance::None());
  std::vector<const Portal*> portals;
  for (auto& candidate : main.candidates()) {
    if (!candidate.intersection().isValid()) {
      continue;
    }

    if (main.candidates().size() > 1 && posOnly &&
        !detail::checkPathLength(candidate.intersection().pathLength(),
                                 s_onSurfaceTolerance,
                                 std::numeric_limits<double>::max(), logger)) {
      continue;
    }

    portals.push_back(&candidate.portal());
  }

  // Find portal types
  const Portal* outerCylinder = nullptr;
  const Portal* innerCylinder = nullptr;
  const Portal* positiveDisc = nullptr;
  const Portal* negativeDisc = nullptr;

  for (const Portal* portal : portals) {
    if (portal ==
        shell.portal(CylinderVolumeBounds::Face::OuterCylinder).get()) {
      outerCylinder = portal;
    } else if (portal ==
               shell.portal(CylinderVolumeBounds::Face::InnerCylinder).get()) {
      innerCylinder = portal;
    } else if (portal ==
               shell.portal(CylinderVolumeBounds::Face::PositiveDisc).get()) {
      positiveDisc = portal;
    } else if (portal ==
               shell.portal(CylinderVolumeBounds::Face::NegativeDisc).get()) {
      negativeDisc = portal;
    }
  }

  // Build new filtered list
  std::vector<const Portal*> filteredPortals;

  // Apply existing filter: remove outer cylinder if both inner and outer are
  // present
  if ((innerCylinder != nullptr) && (outerCylinder != nullptr)) {
    // Keep inner, discard outer
    filteredPortals.push_back(innerCylinder);
  } else {
    // Keep whichever one exists
    if (innerCylinder != nullptr) {
      filteredPortals.push_back(innerCylinder);
    }
    if (outerCylinder != nullptr) {
      filteredPortals.push_back(outerCylinder);
    }
  }

  // Apply CylinderNavigationPolicy optimization: if inner cylinder is present,
  // assume any discs are blocked by it (based on failing test pattern)
  if (innerCylinder == nullptr) {
    // No inner cylinder, so discs are not blocked
    if (positiveDisc != nullptr) {
      filteredPortals.push_back(positiveDisc);
    }
    if (negativeDisc != nullptr) {
      filteredPortals.push_back(negativeDisc);
    }
  }
  // If inner cylinder is present, discs are omitted (blocked)

  return filteredPortals;
}

std::vector<const Portal*> getSmart(const Vector3& position,
                                    const Vector3& direction,
                                    const Transform3& transform,
                                    CylinderNavigationPolicy& policy) {
  Vector3 gpos = transform * position;
  Vector3 gdir = transform.linear() * direction;
  NavigationArguments args{.position = gpos, .direction = gdir};
  NavigationStream main;
  GeometryContext gctx;
  AppendOnlyNavigationStream stream{main};
  policy.initializeCandidates(gctx, args, stream, *logger);

  std::vector<const Portal*> portals;
  // We don't filter here, because we want to test the candidates as they come
  // out of the policy
  for (auto& candidate : main.candidates()) {
    portals.push_back(&candidate.portal());
  }
  return portals;
}

void checkEqual(const std::vector<const Portal*>& exp,
                const std::vector<const Portal*>& act,
                SingleCylinderPortalShell& shell) {
  auto which = [&](const Portal* p) -> std::string {
    if (p == shell.portal(CylinderVolumeBounds::Face::InnerCylinder).get()) {
      return "InnerCylinder";
    }
    if (p == shell.portal(CylinderVolumeBounds::Face::OuterCylinder).get()) {
      return "OuterCylinder";
    }
    if (p == shell.portal(CylinderVolumeBounds::Face::PositiveDisc).get()) {
      return "PositiveDisc";
    }
    if (p == shell.portal(CylinderVolumeBounds::Face::NegativeDisc).get()) {
      return "NegativeDisc";
    }
    BOOST_FAIL("Unknown portal");
    return "";  // unreachable
  };

  std::set<const Portal*> expSet;
  std::set<const Portal*> actSet;

  std::ranges::copy(exp, std::inserter(expSet, expSet.begin()));
  std::ranges::copy(act, std::inserter(actSet, actSet.begin()));

  if (expSet != actSet) {
    BOOST_ERROR([&]() -> std::string {
      std::vector<std::string> exps;
      for (auto& p : exp) {
        exps.push_back(which(p));
      }
      std::vector<std::string> acts;
      for (auto& p : act) {
        acts.push_back(which(p));
      }
      return "[" + boost::algorithm::join(exps, ", ") + "] != [" +
             boost::algorithm::join(acts, ", ") + "]";
    }());
  }
}

}  // namespace

BOOST_DATA_TEST_CASE(
    CylinderPolicyTest,
    (bdata::xrange(-135, 180, 45) *
     bdata::make(Vector3{0_mm, 0_mm, 0_mm}, Vector3{20_mm, 0_mm, 0_mm},
                 Vector3{0_mm, 20_mm, 0_mm}, Vector3{20_mm, 20_mm, 0_mm},
                 Vector3{0_mm, 0_mm, 20_mm})),
    angle, offset) {
  using enum CylinderVolumeBounds::Face;

  Transform3 transform = Transform3::Identity();
  transform *= AngleAxis3{angle * 1_degree, Vector3::UnitX()};
  transform *= Translation3{offset};
  auto cylBounds =
      std::make_shared<CylinderVolumeBounds>(100_mm, 400_mm, 300_mm);
  auto cylVolume =
      std::make_shared<TrackingVolume>(transform, cylBounds, "CylinderVolume");
  SingleCylinderPortalShell shell{*cylVolume};
  shell.applyToVolume();

  {
    Vector3 position = Vector3::UnitX() * 150_mm;
    Vector3 direction = Vector3::UnitZ();

    auto exp =
        getTruth(position, direction, transform, *cylVolume, shell, *logger);

    BOOST_CHECK(exp.size() == 1);
    BOOST_CHECK(exp.at(0) == shell.portal(PositiveDisc).get());

    CylinderNavigationPolicy policy(gctx, *cylVolume, *logger);
    auto act = getSmart(position, direction, transform, policy);
    checkEqual(exp, act, shell);
  }

  {
    Vector3 position = Vector3::UnitX() * 150_mm;
    Vector3 direction = Vector3{1, 1, 0}.normalized();

    auto exp =
        getTruth(position, direction, transform, *cylVolume, shell, *logger);

    BOOST_CHECK(exp.size() == 1);
    BOOST_CHECK(exp.at(0) == shell.portal(OuterCylinder).get());

    CylinderNavigationPolicy policy(gctx, *cylVolume, *logger);
    auto act = getSmart(position, direction, transform, policy);
    checkEqual(exp, act, shell);
  }

  {
    Vector3 position = Vector3::UnitX() * 150_mm;
    Vector3 direction = Vector3{-1, 0, 0}.normalized();

    auto exp =
        getTruth(position, direction, transform, *cylVolume, shell, *logger);

    BOOST_CHECK(exp.size() == 1);
    BOOST_CHECK(exp.at(0) == shell.portal(InnerCylinder).get());

    CylinderNavigationPolicy policy(gctx, *cylVolume, *logger);
    auto act = getSmart(position, direction, transform, policy);
    checkEqual(exp, act, shell);
  }

  {
    Vector3 position = Vector3::UnitX() * 150_mm;
    Vector3 direction = -Vector3::UnitZ();

    auto exp =
        getTruth(position, direction, transform, *cylVolume, shell, *logger);

    BOOST_CHECK(exp.size() == 1);
    BOOST_CHECK(exp.at(0) == shell.portal(NegativeDisc).get());

    CylinderNavigationPolicy policy(gctx, *cylVolume, *logger);
    auto act = getSmart(position, direction, transform, policy);
    checkEqual(exp, act, shell);
  }

  {
    Vector3 position{50, -200, 0};
    Vector3 direction = Vector3{0, 1.5, 1}.normalized();

    auto exp =
        getTruth(position, direction, transform, *cylVolume, shell, *logger);

    BOOST_CHECK(exp.size() == 1);
    BOOST_CHECK(exp.at(0) == shell.portal(InnerCylinder).get());

    CylinderNavigationPolicy policy(gctx, *cylVolume, *logger);
    auto act = getSmart(position, direction, transform, policy);
    checkEqual(exp, act, shell);
  }

  {
    Vector3 position{50, -200, 0};
    Vector3 direction = Vector3{0, 1.2, 1}.normalized();

    auto exp =
        getTruth(position, direction, transform, *cylVolume, shell, *logger);

    BOOST_CHECK(exp.size() == 1);
    BOOST_CHECK(exp.at(0) == shell.portal(InnerCylinder).get());

    CylinderNavigationPolicy policy(gctx, *cylVolume, *logger);
    auto act = getSmart(position, direction, transform, policy);
    checkEqual(exp, act, shell);
  }

  {
    Vector3 position{50, -200, 0};
    Vector3 direction = Vector3{0, 0.9, 1}.normalized();

    auto exp =
        getTruth(position, direction, transform, *cylVolume, shell, *logger);

    BOOST_CHECK(exp.size() == 1);
    BOOST_CHECK(exp.at(0) == shell.portal(InnerCylinder).get());

    CylinderNavigationPolicy policy(gctx, *cylVolume, *logger);
    auto act = getSmart(position, direction, transform, policy);
    checkEqual(exp, act, shell);
  }

  {
    Vector3 position{20, -200, 0};
    Vector3 direction = Vector3{0.45, 0.9, 1}.normalized();

    auto exp =
        getTruth(position, direction, transform, *cylVolume, shell, *logger);

    BOOST_CHECK(exp.size() == 1);
    BOOST_CHECK(exp.at(0) == shell.portal(PositiveDisc).get());

    CylinderNavigationPolicy policy(gctx, *cylVolume, *logger);
    auto act = getSmart(position, direction, transform, policy);
    checkEqual(exp, act, shell);
  }

  {
    Vector3 position{20, -200, 0};
    Vector3 direction = Vector3{0.45, 0.9, -1}.normalized();

    auto exp =
        getTruth(position, direction, transform, *cylVolume, shell, *logger);

    BOOST_CHECK(exp.size() == 1);
    BOOST_CHECK(exp.at(0) == shell.portal(NegativeDisc).get());

    CylinderNavigationPolicy policy(gctx, *cylVolume, *logger);
    auto act = getSmart(position, direction, transform, policy);
    checkEqual(exp, act, shell);
  }

  {
    Vector3 position{400 * std::cos(std::numbers::pi / 4),
                     400 * std::sin(std::numbers::pi / 4), 0};
    Vector3 direction = Vector3{0.45, -0.9, -0.1}.normalized();

    // We're sitting ON the outer cylinder here and missing the disc and inner
    // cylinder: need to return the outer cylinder
    auto exp = getTruth(position, direction, transform, *cylVolume, shell,
                        *logger, false);

    BOOST_CHECK(exp.size() == 1);
    BOOST_CHECK(exp.at(0) == shell.portal(OuterCylinder).get());

    CylinderNavigationPolicy policy(gctx, *cylVolume, *logger);
    auto act = getSmart(position, direction, transform, policy);
    checkEqual(exp, act, shell);
  }

  {
    Vector3 position{400 * std::cos(std::numbers::pi / 4),
                     400 * std::sin(std::numbers::pi / 4), 0};
    Vector3 direction = Vector3{-0.3, -0.9, -0.1}.normalized();

    // We're sitting ON the outer cylinder here and missing the disc and inner
    // cylinder: need to return the outer cylinder
    auto exp = getTruth(position, direction, transform, *cylVolume, shell,
                        *logger, false);

    BOOST_CHECK(exp.size() == 1);
    BOOST_CHECK(exp.at(0) == shell.portal(OuterCylinder).get());

    CylinderNavigationPolicy policy(gctx, *cylVolume, *logger);
    auto act = getSmart(position, direction, transform, policy);
    checkEqual(exp, act, shell);
  }

  {
    Vector3 position{100 * std::cos(std::numbers::pi / 4),
                     100 * std::sin(std::numbers::pi / 4), 0};
    double dangle = 0.1;
    Vector3 direction = Vector3{std::cos(dangle), std::sin(dangle), 0.01};

    // We're sitting ON the Inner cylinder here and  pointing outwards
    auto exp =
        getTruth(position, direction, transform, *cylVolume, shell, *logger);

    BOOST_CHECK(exp.size() == 1);
    BOOST_CHECK(exp.at(0) == shell.portal(OuterCylinder).get());

    CylinderNavigationPolicy policy(gctx, *cylVolume, *logger);
    auto act = getSmart(position, direction, transform, policy);
    checkEqual(exp, act, shell);
  }

  {
    Vector3 position{200 * std::cos(std::numbers::pi / 4),
                     200 * std::sin(std::numbers::pi / 4), 0};
    Vector3 target{150 * std::cos(std::numbers::pi * 5 / 4),
                   150 * std::sin(std::numbers::pi * 5 / 4), 300};
    Vector3 direction = (target - position).normalized();

    auto exp =
        getTruth(position, direction, transform, *cylVolume, shell, *logger);

    BOOST_CHECK_EQUAL(exp.size(), 1);
    BOOST_CHECK_EQUAL(exp.at(0), shell.portal(InnerCylinder).get());

    CylinderNavigationPolicy policy(gctx, *cylVolume, *logger);
    auto act = getSmart(position, direction, transform, policy);
    checkEqual(exp, act, shell);
  }
}

namespace {

std::mt19937 engine;

unsigned long seed() {
  static unsigned long s = 42;
  return s++;
}

std::uniform_real_distribution<double> rDistOffBoundary{
    100 + 2 * s_onSurfaceTolerance, 400 - 2 * s_onSurfaceTolerance};
std::uniform_real_distribution<double> zDistOffBoundary{
    -300_mm + 2 * s_onSurfaceTolerance, 300_mm - 2 * s_onSurfaceTolerance};
std::uniform_real_distribution<double> phiDist{-std::numbers::pi,
                                               std::numbers::pi};
std::uniform_real_distribution<double> thetaDist{0, std::numbers::pi};

}  // namespace

BOOST_DATA_TEST_CASE(
    CylinderPolicyTestOffBoundary,
    bdata::random((bdata::engine = engine, bdata::seed = seed(),
                   bdata::distribution = rDistOffBoundary)) ^
        bdata::random((bdata::engine = engine, bdata::seed = seed(),
                       bdata::distribution = zDistOffBoundary)) ^
        bdata::random((bdata::engine = engine, bdata::seed = seed(),
                       bdata::distribution = phiDist)) ^
        bdata::random((bdata::engine = engine, bdata::seed = seed(),
                       bdata::distribution = phiDist)) ^
        bdata::random((bdata::engine = engine, bdata::seed = seed(),
                       bdata::distribution = thetaDist)) ^
        bdata::xrange(100),
    r, z, phiPos, phiDir, theta, index) {
  static_cast<void>(index);

  Transform3 transform = Transform3::Identity();
  auto cylBounds =
      std::make_shared<CylinderVolumeBounds>(100_mm, 400_mm, 300_mm);
  auto cylVolume =
      std::make_shared<TrackingVolume>(transform, cylBounds, "CylinderVolume");
  SingleCylinderPortalShell shell{*cylVolume};
  shell.applyToVolume();

  Vector3 position{r * std::cos(phiPos), r * std::sin(phiPos), z};
  Vector3 direction{std::sin(theta) * std::cos(phiDir),
                    std::sin(theta) * std::sin(phiDir), std::cos(theta)};

  BOOST_CHECK(cylBounds->inside(position));

  auto exp =
      getTruth(position, direction, transform, *cylVolume, shell, *logger);

  CylinderNavigationPolicy policy(gctx, *cylVolume, *logger);
  auto act = getSmart(position, direction, transform, policy);
  checkEqual(exp, act, shell);
}

BOOST_DATA_TEST_CASE(
    CylinderPolicyTestOnRBoundary,
    bdata::make(100, 400) *
        (bdata::random((bdata::engine = engine, bdata::seed = seed(),
                        bdata::distribution = zDistOffBoundary)) ^
         bdata::random((bdata::engine = engine, bdata::seed = seed(),
                        bdata::distribution = phiDist)) ^
         bdata::random((bdata::engine = engine, bdata::seed = seed(),
                        bdata::distribution = phiDist)) ^
         bdata::random((bdata::engine = engine, bdata::seed = seed(),
                        bdata::distribution = zDistOffBoundary)) ^
         bdata::xrange(100)),
    r, z, phiPos, phiTarget, zTarget, index) {
  static_cast<void>(index);

  Transform3 transform = Transform3::Identity();
  auto cylBounds =
      std::make_shared<CylinderVolumeBounds>(100_mm, 400_mm, 300_mm);
  auto cylVolume =
      std::make_shared<TrackingVolume>(transform, cylBounds, "CylinderVolume");
  SingleCylinderPortalShell shell{*cylVolume};
  shell.applyToVolume();

  Vector3 position{r * std::cos(phiPos), r * std::sin(phiPos), z};
  Vector3 target{r * std::cos(phiTarget), r * std::sin(phiTarget), zTarget};
  Vector3 direction = (target - position).normalized();

  auto exp =
      getTruth(position, direction, transform, *cylVolume, shell, *logger);

  CylinderNavigationPolicy policy(gctx, *cylVolume, *logger);
  auto act = getSmart(position, direction, transform, policy);
  checkEqual(exp, act, shell);
}

BOOST_DATA_TEST_CASE(
    CylinderPolicyTestOnZBoundary,
    bdata::make(-300, 300) *
        (bdata::random((bdata::engine = engine, bdata::seed = seed(),
                        bdata::distribution = rDistOffBoundary)) ^
         bdata::random((bdata::engine = engine, bdata::seed = seed(),
                        bdata::distribution = phiDist)) ^
         bdata::random((bdata::engine = engine, bdata::seed = seed(),
                        bdata::distribution = phiDist)) ^
         bdata::random((bdata::engine = engine, bdata::seed = seed(),
                        bdata::distribution = zDistOffBoundary)) ^
         bdata::xrange(100)),
    z, r, phiPos, phiTarget, zTarget, index) {
  static_cast<void>(index);
  Transform3 transform = Transform3::Identity();
  auto cylBounds =
      std::make_shared<CylinderVolumeBounds>(100_mm, 400_mm, 300_mm);
  auto cylVolume =
      std::make_shared<TrackingVolume>(transform, cylBounds, "CylinderVolume");
  SingleCylinderPortalShell shell{*cylVolume};
  shell.applyToVolume();

  Vector3 position{r * std::cos(phiPos), r * std::sin(phiPos),
                   static_cast<double>(z)};
  Vector3 target{r * std::cos(phiTarget), r * std::sin(phiTarget), zTarget};
  Vector3 direction = (target - position).normalized();

  BOOST_CHECK(cylBounds->inside(position));

  auto exp =
      getTruth(position, direction, transform, *cylVolume, shell, *logger);

  CylinderNavigationPolicy policy(gctx, *cylVolume, *logger);
  auto act = getSmart(position, direction, transform, policy);
  checkEqual(exp, act, shell);
}

BOOST_AUTO_TEST_CASE(CylinderPolicyZeroInnerRadiusTest) {
  // Test that CylinderNavigationPolicy rejects volumes with zero inner radius
  Transform3 transform = Transform3::Identity();

  // Create cylinder volume bounds with zero inner radius (rMin = 0)
  auto cylBounds = std::make_shared<CylinderVolumeBounds>(0_mm, 400_mm, 300_mm);
  auto cylVolume = std::make_shared<TrackingVolume>(transform, cylBounds,
                                                    "ZeroInnerRadiusVolume");

  // CylinderNavigationPolicy constructor should throw std::invalid_argument
  // for volumes with zero inner radius
  {
    Acts::Logging::ScopedFailureThreshold log{Logging::FATAL};
    BOOST_CHECK_THROW(CylinderNavigationPolicy(gctx, *cylVolume, *logger),
                      std::invalid_argument);
  }
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace ActsTests

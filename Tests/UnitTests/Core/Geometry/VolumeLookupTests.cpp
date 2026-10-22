// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Geometry/Blueprint.hpp"
#include "Acts/Geometry/CuboidVolumeBounds.hpp"
#include "Acts/Geometry/CylinderVolumeBounds.hpp"
#include "Acts/Geometry/LinearVolumeLookup.hpp"
#include "Acts/Geometry/RzGridVolumeLookup.hpp"
#include "Acts/Geometry/StaticBlueprintNode.hpp"
#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/Geometry/TrackingVolume.hpp"
#include "Acts/Geometry/VolumePlacementBase.hpp"

#include <limits>
#include <random>
#include <sstream>
#include <type_traits>

using namespace Acts;

static_assert(!std::is_copy_constructible_v<TrackingGeometry>);
static_assert(!std::is_copy_assignable_v<TrackingGeometry>);
namespace {
auto gctx = GeometryContext::dangerouslyDefaultConstruct();
std::shared_ptr<TrackingVolume> makeWorld() {
  auto world = std::make_shared<TrackingVolume>(
      Transform3::Identity(),
      std::make_shared<CylinderVolumeBounds>(0., 100., 100.), "world");
  // Overlapping siblings and gaps exercise priority and parent fallback.
  for (int i = 0; i < 16; ++i) {
    auto child = std::make_unique<TrackingVolume>(
        Transform3::Identity(),
        std::make_shared<CylinderVolumeBounds>(i * 4., i * 4. + 6., 20. + i),
        "shell");
    if (i == 2) {
      child->addVolume(std::make_unique<TrackingVolume>(
          Transform3::Identity(),
          std::make_shared<CylinderVolumeBounds>(9., 10., 5.), "nested"));
    }
    world->addVolume(std::move(child));
  }
  unsigned id = 0;
  world->apply([&](TrackingVolume& v) {
    v.assignGeometryId(GeometryIdentifier{}.withVolume(++id));
  });
  return world;
}
class MovingPlacement : public VolumePlacementBase {
 public:
  Transform3 transform = Transform3::Identity();
  Transform3 inverse = Transform3::Identity();
  const Transform3& localToGlobalTransform(
      const GeometryContext& /*gctx*/) const override {
    return transform;
  }
  const Transform3& portalLocalToGlobal(
      const GeometryContext& /*gctx*/,
      std::size_t /*portalIdx*/) const override {
    return transform;
  }
  const Transform3& globalToLocalTransform(
      const GeometryContext& /*gctx*/) const override {
    return inverse;
  }
};
}  // namespace

BOOST_AUTO_TEST_CASE(RzLookupMatchesLinear) {
  auto world = makeWorld();
  const LinearVolumeLookup linear(*world, {});
  TrackingGeometry geometry(world, nullptr, {}, getDummyLogger(), false);
  std::mt19937 rng(42);
  std::uniform_real_distribution<double> uniform(-110., 110.);
  {
    RzGridVolumeLookup::Config cfg;
    geometry.highestTrackingVolume()->setVolumeLookup(
        std::make_unique<RzGridVolumeLookup>(gctx, *world, cfg));
    for (double tolerance : {-1e-4, 0., 1e-6, cfg.maxTolerance, 0.5}) {
      for (int i = 0; i < 3000; ++i) {
        Vector3 p(uniform(rng), uniform(rng), uniform(rng));
        BOOST_CHECK(*geometry.resolveLowestTrackingVolume(gctx, p, std::nullopt,
                                                          nullptr, tolerance) ==
                    linear.lowestTrackingVolume(gctx, p, tolerance));
      }
      for (int i = 0; i < 16; ++i) {
        for (double delta : {-2e-3, -1e-3, -1e-6, 0., 1e-6, 1e-3, 2e-3}) {
          for (Vector3 p : {Vector3(i * 4. + delta, 0., 0.),
                            Vector3(i * 4. + 6. + delta, 0., 0.),
                            Vector3(i * 4. + 2., 0., 20. + i + delta)}) {
            BOOST_CHECK(*geometry.resolveLowestTrackingVolume(
                            gctx, p, std::nullopt, nullptr, tolerance) ==
                        linear.lowestTrackingVolume(gctx, p, tolerance));
          }
        }
      }
    }
    // First sibling wins even when a later sibling has a deeper matching child.
    BOOST_CHECK(*geometry.resolveLowestTrackingVolume(gctx, {9.5, 0., 0.}) ==
                &world->volumes().at(1));
    geometry.highestTrackingVolume()->setVolumeLookup(nullptr);
    BOOST_CHECK(*geometry.resolveLowestTrackingVolume(gctx, {9.5, 0., 0.}) ==
                &world->volumes().at(1));
    BOOST_CHECK(geometry.highestTrackingVolume()->volumeLookup() == nullptr);
  }
}

BOOST_AUTO_TEST_CASE(RzLookupFallbacks) {
  {
    auto world = makeWorld();
    RzGridVolumeLookup::Config cfg;
    {
      cfg.maxGridCells = 1;  // Force bounded-memory fallback.
    }
    RzGridVolumeLookup limited(gctx, *world, cfg);
    BOOST_CHECK(limited.lowestTrackingVolume(gctx, {40., 0., 0.}, 0.) ==
                world->lowestTrackingVolume(gctx, {40., 0., 0.}, 0.));
    world->addVolume(std::make_unique<TrackingVolume>(
        Transform3::Identity(),
        std::make_shared<CuboidVolumeBounds>(80., 80., 80.), "box"));
    RzGridVolumeLookup mixed(gctx, *world, cfg);
    BOOST_CHECK(mixed.lowestTrackingVolume(gctx, {70., 0., 0.}, 0.) ==
                world->lowestTrackingVolume(gctx, {70., 0., 0.}, 0.));

    MovingPlacement placement;
    auto movingWorld = std::make_shared<TrackingVolume>(
        Transform3::Identity(),
        std::make_shared<CylinderVolumeBounds>(0., 100., 100.));
    movingWorld->addVolume(std::make_unique<TrackingVolume>(
        placement, std::make_shared<CylinderVolumeBounds>(0., 5., 5.)));
    RzGridVolumeLookup moving(gctx, *movingWorld, cfg);
    placement.transform.translation().x() = 50.;
    placement.inverse = placement.transform.inverse();
    BOOST_CHECK(moving.lowestTrackingVolume(gctx, {50., 0., 0.}, 0.) ==
                &movingWorld->volumes().at(0));
  }
}

BOOST_AUTO_TEST_CASE(BlueprintSelectsVolumeLookup) {
  for (int type : {0, 1}) {
    Blueprint blueprint{Blueprint::Config{}};
    blueprint.addChild(
        std::make_shared<StaticBlueprintNode>(std::make_unique<TrackingVolume>(
            Transform3::Identity(),
            std::make_shared<CylinderVolumeBounds>(0., 10., 10.), "cylinder")));
    BlueprintOptions options;
    blueprint.setVolumeLookupFactory(
        [type](const GeometryContext& context, const TrackingVolume& volume,
               const Logger& logger) -> std::unique_ptr<IVolumeLookup> {
          if (type == 0) {
            return std::make_unique<LinearVolumeLookup>(
                volume, LinearVolumeLookup::Config{}, logger);
          }
          return std::make_unique<RzGridVolumeLookup>(
              context, volume, RzGridVolumeLookup::Config{}, logger);
        });
    auto geometry = blueprint.construct(options, gctx);
    geometry->highestTrackingVolume()->apply([&](TrackingVolume& volume) {
      if (&volume != geometry->highestTrackingVolume()) {
        BOOST_CHECK(volume.volumeLookup() == nullptr);
        return;
      }
      BOOST_REQUIRE(volume.volumeLookup() != nullptr);
      BOOST_CHECK((dynamic_cast<const LinearVolumeLookup*>(
                       volume.volumeLookup()) != nullptr) == (type == 0));
      BOOST_CHECK((dynamic_cast<const RzGridVolumeLookup*>(
                       volume.volumeLookup()) != nullptr) == (type == 1));
    });
    BOOST_CHECK(*geometry->resolveLowestTrackingVolume(gctx, {5., 0., 0.}) !=
                nullptr);
  }
}

BOOST_AUTO_TEST_CASE(InvalidVolumeLookupConfig) {
  auto world = makeWorld();
  RzGridVolumeLookup::Config cfg;
  cfg.maxTolerance = -1;
  BOOST_CHECK_THROW(RzGridVolumeLookup(gctx, *world, cfg),
                    std::invalid_argument);
  cfg.maxTolerance = 0;
  cfg.maxGridCells = 0;
  BOOST_CHECK_THROW(RzGridVolumeLookup(gctx, *world, cfg),
                    std::invalid_argument);
}

BOOST_AUTO_TEST_CASE(InvalidGridConfig) {
  auto world = makeWorld();
  RzGridVolumeLookup::Config cfg;
  cfg.maxGridCells = 0;
  BOOST_CHECK_THROW(RzGridVolumeLookup(gctx, *world, cfg),
                    std::invalid_argument);
}

BOOST_AUTO_TEST_CASE(LinearLookupWarnings) {
  // These assertions deliberately exercise warning diagnostics.
#ifdef ACTS_ENABLE_LOG_FAILURE_THRESHOLD
  Logging::ScopedFailureThreshold threshold{Logging::ERROR};
#endif
  auto world = makeWorld();
  std::ostringstream output;
  auto logger = getDefaultLogger("LookupTest", Logging::WARNING, &output);
  LinearVolumeLookup atThreshold(*world, LinearVolumeLookup::Config{16},
                                 *logger);
  BOOST_CHECK(output.str().empty());
  LinearVolumeLookup aboveThreshold(*world, LinearVolumeLookup::Config{15},
                                    *logger);
  BOOST_CHECK(output.str().find("scans up to 16 siblings") !=
              std::string::npos);
  output.str("");
  LinearVolumeLookup suppressed(*world, LinearVolumeLookup::Config{0}, *logger);
  BOOST_CHECK(output.str().empty());
  for (int i = 16; i < 65; ++i) {
    world->addVolume(std::make_unique<TrackingVolume>(
        Transform3::Identity(),
        std::make_shared<CylinderVolumeBounds>(0., 5., 5.), "extra"));
  }
  unsigned id = 0;
  world->apply([&](TrackingVolume& v) {
    v.assignGeometryId(GeometryIdentifier{}.withVolume(++id));
  });
  // A configured accelerator does not emit a transient linear warning.
  auto lookup = std::make_unique<RzGridVolumeLookup>(
      gctx, *world, RzGridVolumeLookup::Config{}, *logger);
  world->setVolumeLookup(std::move(lookup));
  TrackingGeometry accelerated(world, nullptr, {}, *logger, false);
  BOOST_CHECK(output.str().empty());
  world->setVolumeLookup(nullptr);
  BOOST_CHECK(output.str().empty());
  BOOST_CHECK(world->volumeLookup() == nullptr);
  world->setVolumeLookup(std::make_unique<LinearVolumeLookup>(
      *world, LinearVolumeLookup::Config{}, *logger));
  BOOST_CHECK(output.str().find("scans up to 65 siblings") !=
              std::string::npos);
}

BOOST_AUTO_TEST_CASE(AcceleratorFallbackWarnings) {
  // These assertions deliberately exercise warning diagnostics.
#ifdef ACTS_ENABLE_LOG_FAILURE_THRESHOLD
  Logging::ScopedFailureThreshold threshold{Logging::ERROR};
#endif
  {
    std::ostringstream output;
    auto logger = getDefaultLogger("LookupTest", Logging::WARNING, &output);
    auto world = std::make_shared<TrackingVolume>(
        Transform3::Identity(),
        std::make_shared<CylinderVolumeBounds>(0., 100., 100.));
    auto container = std::make_unique<TrackingVolume>(
        Transform3::Identity(),
        std::make_shared<CylinderVolumeBounds>(0., 50., 50.));
    container->addVolume(std::make_unique<TrackingVolume>(
        Transform3::Identity(),
        std::make_shared<CuboidVolumeBounds>(2., 2., 2.)));
    world->addVolume(std::move(container));
    world->addVolume(std::make_unique<TrackingVolume>(
        Transform3::Identity(),
        std::make_shared<CylinderVolumeBounds>(50., 90., 50.)));
    RzGridVolumeLookup::Config cfg;
    cfg.fallbackWarningFraction = 0.3;
    RzGridVolumeLookup warns(gctx, *world, cfg, *logger);
    BOOST_CHECK(output.str().find("1/3 volumes in the indexed scope") !=
                std::string::npos);
    BOOST_CHECK(output.str().find("unsupported shape/orientation=1") !=
                std::string::npos);
    output.str("");
    cfg.fallbackWarningFraction = 1. / 3.;
    RzGridVolumeLookup equal(gctx, *world, cfg, *logger);
    BOOST_CHECK(output.str().empty());
    cfg.fallbackWarningFraction = 1.;
    RzGridVolumeLookup suppressed(gctx, world->volumes().at(0), cfg, *logger);
    BOOST_CHECK(output.str().empty());

    MovingPlacement placement;
    auto movingWorld = std::make_shared<TrackingVolume>(
        Transform3::Identity(),
        std::make_shared<CylinderVolumeBounds>(0., 100., 100.));
    movingWorld->addVolume(std::make_unique<TrackingVolume>(
        placement, std::make_shared<CylinderVolumeBounds>(0., 10., 10.)));
    cfg.fallbackWarningFraction = 0.;
    RzGridVolumeLookup moving(gctx, *movingWorld, cfg, *logger);
    BOOST_CHECK(output.str().find("1/1 volumes in the indexed scope") !=
                std::string::npos);
    BOOST_CHECK(output.str().find("alignable=1") != std::string::npos);
    output.str("");
    {
      cfg.maxGridCells = 1;
      auto largeWorld = makeWorld();
      RzGridVolumeLookup limited(gctx, *largeWorld, cfg, *logger);
      BOOST_CHECK(output.str().find("grid cell limit=1") != std::string::npos);
    }
    for (double invalid :
         {-0.1, 1.1, std::numeric_limits<double>::quiet_NaN()}) {
      cfg.fallbackWarningFraction = invalid;
      BOOST_CHECK_THROW(RzGridVolumeLookup(gctx, *world, cfg, *logger),
                        std::invalid_argument);
    }
  }
}

BOOST_AUTO_TEST_CASE(BlueprintLookupWarningSettings) {
  // These assertions deliberately exercise warning diagnostics.
#ifdef ACTS_ENABLE_LOG_FAILURE_THRESHOLD
  Logging::ScopedFailureThreshold threshold{Logging::ERROR};
#endif
  std::ostringstream output;
  auto logger = getDefaultLogger("BlueprintTest", Logging::WARNING, &output);
  for (int type : {0, 1}) {
    Blueprint blueprint{Blueprint::Config{}};
    blueprint.addChild(
        std::make_shared<StaticBlueprintNode>(std::make_unique<TrackingVolume>(
            Transform3::Identity(),
            std::make_shared<CuboidVolumeBounds>(10., 10., 10.), "cube")));
    BlueprintOptions options;
    blueprint.setVolumeLookupFactory(
        [type](const GeometryContext& context, const TrackingVolume& volume,
               const Logger& factoryLogger) -> std::unique_ptr<IVolumeLookup> {
          if (type == 0) {
            return std::make_unique<LinearVolumeLookup>(
                volume, LinearVolumeLookup::Config{0}, factoryLogger);
          }
          RzGridVolumeLookup::Config config;
          config.fallbackWarningFraction = 1;
          return std::make_unique<RzGridVolumeLookup>(context, volume, config,
                                                      factoryLogger);
        });
    auto geometry = blueprint.construct(options, gctx, *logger);
    BOOST_CHECK(output.str().find("volume lookup") == std::string::npos);
    if (auto* linear = dynamic_cast<const LinearVolumeLookup*>(
            geometry->highestTrackingVolume()->volumeLookup())) {
      BOOST_CHECK_EQUAL(linear->config().warningThreshold, 0);
    }
    if (auto* grid = dynamic_cast<const RzGridVolumeLookup*>(
            geometry->highestTrackingVolume()->volumeLookup())) {
      BOOST_CHECK_EQUAL(grid->config().fallbackWarningFraction, 1);
    }
  }
}

BOOST_AUTO_TEST_CASE(RecursiveLookupUsesEachVolumesPolicy) {
  class CountingLookup final : public IVolumeLookup {
   public:
    CountingLookup(const TrackingVolume& root, unsigned& calls)
        : m_calls(calls), m_linear(root, {}) {}
    const TrackingVolume* lowestTrackingVolume(
        const GeometryContext& context, const Vector3& position,
        double tolerance) const override {
      ++m_calls;
      return m_linear.lowestTrackingVolume(context, position, tolerance);
    }

   private:
    unsigned& m_calls;
    LinearVolumeLookup m_linear;
  };
  {
    auto world = makeWorld();
    unsigned calls = 0;
    auto& child = world->volumes().at(1);
    child.setVolumeLookup(std::make_unique<CountingLookup>(child, calls));
    TrackingGeometry geometry(world, nullptr, {}, getDummyLogger(), false);
    BOOST_CHECK(*geometry.resolveLowestTrackingVolume(gctx, {9.5, 0., 0.}) ==
                &child);
    BOOST_CHECK_EQUAL(calls, 1);
    BOOST_CHECK(world->lowestTrackingVolume(gctx, {9.5, 0., 0.}) == &child);
    BOOST_CHECK_EQUAL(calls, 2);
    world->setVolumeLookup(std::make_unique<RzGridVolumeLookup>(
        gctx, *world, RzGridVolumeLookup::Config{}));
    BOOST_CHECK(world->lowestTrackingVolume(gctx, {9.5, 0., 0.}) == &child);
    BOOST_CHECK_EQUAL(calls, 3);
    // An unsupported sibling is checked linearly alongside indexed candidates.
    // The chosen child's policy must still run. Reset only the parent, leaving
    // the child intact.
    world->addVolume(std::make_unique<TrackingVolume>(
        Transform3::Identity(),
        std::make_shared<CuboidVolumeBounds>(1., 1., 1.)));
    world->setVolumeLookup(std::make_unique<RzGridVolumeLookup>(
        gctx, *world, RzGridVolumeLookup::Config{}, getDummyLogger()));
    BOOST_CHECK(*geometry.resolveLowestTrackingVolume(gctx, {9.5, 0., 0.}) ==
                &child);
    BOOST_CHECK_EQUAL(calls, 4);
    world->setVolumeLookup(nullptr);
    BOOST_CHECK(world->lowestTrackingVolume(gctx, {9.5, 0., 0.}) == &child);
    BOOST_CHECK_EQUAL(calls, 5);
    child.setVolumeLookup(nullptr);
    BOOST_CHECK(child.volumeLookup() == nullptr);
  }
}

BOOST_AUTO_TEST_CASE(LookupBindsToFinalVolumeAddress) {
  auto world = makeWorld();
  TrackingVolume moved(std::move(*world));
  moved.setVolumeLookup(std::make_unique<RzGridVolumeLookup>(
      gctx, moved, RzGridVolumeLookup::Config{}));
  BOOST_CHECK(moved.lowestTrackingVolume(gctx, {90., 0., 0.}) == &moved);
  BOOST_CHECK(moved.lowestTrackingVolume(gctx, {9.5, 0., 0.}) ==
              &moved.volumes().at(1));
  moved.setVolumeLookup(std::make_unique<LinearVolumeLookup>(
      moved, LinearVolumeLookup::Config{}));
  BOOST_CHECK(moved.lowestTrackingVolume(gctx, {90., 0., 0.}) == &moved);
  BOOST_CHECK(moved.lowestTrackingVolume(gctx, {9.5, 0., 0.}) ==
              &moved.volumes().at(1));
}

BOOST_AUTO_TEST_CASE(BlueprintLocalLookupFactories) {
  Blueprint::Config config;
  config.envelope[AxisDirection::AxisR] = {1., 1.};
  config.envelope[AxisDirection::AxisZ] = {1., 1.};
  Blueprint blueprint{config};
  auto region =
      std::make_shared<StaticBlueprintNode>(std::make_unique<TrackingVolume>(
          Transform3::Identity(),
          std::make_shared<CylinderVolumeBounds>(0., 20., 20.), "region"));
  region->addChild(
      std::make_shared<StaticBlueprintNode>(std::make_unique<TrackingVolume>(
          Transform3::Identity(),
          std::make_shared<CylinderVolumeBounds>(0., 10., 10.), "child")));
  blueprint.addChild(region);
  std::vector<std::string> calls;
  auto logger = getDefaultLogger("FactoryTest", Logging::WARNING);
  region->setVolumeLookupFactory([&](const GeometryContext& context,
                                     const TrackingVolume& volume,
                                     const Logger& passedLogger) {
    BOOST_CHECK(&context == &gctx);
    BOOST_CHECK(&passedLogger == logger.get());
    BOOST_CHECK_EQUAL(volume.volumes().size(), 1);
    BOOST_CHECK(volume.volumes().at(0).volumeLookup() == nullptr);
    calls.push_back(volume.volumeName());
    return std::make_unique<LinearVolumeLookup>(
        volume, LinearVolumeLookup::Config{}, passedLogger);
  });
  blueprint.setVolumeLookupFactory([&](const GeometryContext& context,
                                       const TrackingVolume& volume,
                                       const Logger& passedLogger) {
    BOOST_CHECK(&context == &gctx);
    BOOST_CHECK(&passedLogger == logger.get());
    BOOST_CHECK(dynamic_cast<const LinearVolumeLookup*>(
                    volume.volumes().at(0).volumeLookup()) != nullptr);
    calls.push_back(volume.volumeName());
    return std::make_unique<RzGridVolumeLookup>(
        context, volume, RzGridVolumeLookup::Config{}, passedLogger);
  });
  auto geometry = blueprint.construct(BlueprintOptions{}, gctx, *logger);
  BOOST_REQUIRE_EQUAL(calls.size(), 2);
  BOOST_CHECK_EQUAL(calls[0], "region");
  BOOST_CHECK_EQUAL(calls[1], geometry->highestTrackingVolume()->volumeName());
  const auto* result =
      *geometry->resolveLowestTrackingVolume(gctx, {5., 0., 0.});
  BOOST_REQUIRE(result != nullptr);
  BOOST_CHECK_EQUAL(result->volumeName(), "child");
}

BOOST_AUTO_TEST_CASE(BlueprintRejectsNullLookupFactoryResult) {
  for (bool onStatic : {false, true}) {
    Blueprint blueprint{Blueprint::Config{}};
    auto child =
        std::make_shared<StaticBlueprintNode>(std::make_unique<TrackingVolume>(
            Transform3::Identity(),
            std::make_shared<CylinderVolumeBounds>(0., 10., 10.), "child"));
    blueprint.addChild(child);
    VolumeLookupFactory factory = [](const GeometryContext&,
                                     const TrackingVolume&, const Logger&) {
      return std::unique_ptr<IVolumeLookup>{};
    };
    if (onStatic) {
      child->setVolumeLookupFactory(factory);
    } else {
      blueprint.setVolumeLookupFactory(factory);
    }
    BOOST_CHECK_THROW(blueprint.construct(BlueprintOptions{}, gctx),
                      std::invalid_argument);
  }
}

BOOST_AUTO_TEST_CASE(BlueprintDefaultLookupIsImplicit) {
  Blueprint blueprint{Blueprint::Config{}};
  blueprint.addChild(
      std::make_shared<StaticBlueprintNode>(std::make_unique<TrackingVolume>(
          Transform3::Identity(),
          std::make_shared<CylinderVolumeBounds>(0., 10., 10.), "child")));
  auto geometry = blueprint.construct(BlueprintOptions{}, gctx);
  geometry->highestTrackingVolume()->apply([](TrackingVolume& volume) {
    BOOST_CHECK(volume.volumeLookup() == nullptr);
  });
  BOOST_CHECK(*geometry->resolveLowestTrackingVolume(gctx, {5., 0., 0.}) !=
              nullptr);
}

BOOST_AUTO_TEST_CASE(SubtreeLookupPreservesSiblingPriority) {
  auto makeVolume = [](double lo, double hi, const std::string& name) {
    return std::make_unique<TrackingVolume>(
        Transform3::Identity(),
        std::make_shared<CylinderVolumeBounds>(lo, hi, 50.), name);
  };
  auto world = makeVolume(0., 100., "world");
  auto first = makeVolume(0., 20., "first");
  auto inner = makeVolume(5., 15., "inner");
  inner->addVolume(makeVolume(5., 10., "leaf"));
  first->addVolume(std::move(inner));
  world->addVolume(std::move(first));
  auto second = makeVolume(0., 50., "second");
  second->addVolume(makeVolume(10., 14., "later-deeper"));
  world->addVolume(std::move(second));
  const std::vector<std::pair<double, std::string>> points{{7., "leaf"},
                                                           {12., "inner"},
                                                           {18., "first"},
                                                           {30., "second"},
                                                           {80., "world"}};
  {
    world->setVolumeLookup(std::make_unique<RzGridVolumeLookup>(
        gctx, *world, RzGridVolumeLookup::Config{}));
    world->visitVolumes([&](const TrackingVolume* volume) {
      if (volume != world.get()) {
        BOOST_CHECK(volume->volumeLookup() == nullptr);
      }
    });
    for (const auto& [radius, name] : points) {
      const auto* found = world->lowestTrackingVolume(gctx, {radius, 0., 0.});
      BOOST_REQUIRE(found != nullptr);
      BOOST_CHECK_EQUAL(found->volumeName(), name);
    }
    BOOST_CHECK(world->lowestTrackingVolume(gctx, {120., 0., 0.}) == nullptr);
  }
}

BOOST_AUTO_TEST_CASE(SubtreeLookupDelegatesThroughLinearToAnotherIndex) {
  class CountingLinear final : public IVolumeLookup {
   public:
    CountingLinear(const TrackingVolume& root, unsigned& count)
        : m_count(count), m_linear(root, {}) {}
    const TrackingVolume* lowestTrackingVolume(
        const GeometryContext& context, const Vector3& position,
        double tolerance) const override {
      ++m_count;
      return m_linear.lowestTrackingVolume(context, position, tolerance);
    }

   private:
    unsigned& m_count;
    LinearVolumeLookup m_linear;
  };
  auto makeVolume = [](double lo, double hi) {
    return std::make_unique<TrackingVolume>(
        Transform3::Identity(),
        std::make_shared<CylinderVolumeBounds>(lo, hi, 50.));
  };
  auto world = makeVolume(0., 100.);
  auto region = makeVolume(0., 50.);
  auto nested = makeVolume(0., 30.);
  auto leaf = makeVolume(10., 20.);
  const auto* expected = leaf.get();
  nested->addVolume(std::move(leaf));
  nested->setVolumeLookup(std::make_unique<RzGridVolumeLookup>(
      gctx, *nested, RzGridVolumeLookup::Config{}));
  region->addVolume(std::move(nested));
  unsigned calls = 0;
  region->setVolumeLookup(std::make_unique<CountingLinear>(*region, calls));
  world->addVolume(std::move(region));
  {
    std::ostringstream output;
    auto logger = getDefaultLogger("DelegationTest", Logging::WARNING, &output);
    RzGridVolumeLookup::Config config;
    {
      // One boundary entry fits in one regular cell (nine with reserved
      // under-/overflow slots). Indexing through it would exceed
      // this cap and warn, even if fallback happened to return the same result.
      config.maxGridCells = 9;
    }
    world->setVolumeLookup(
        std::make_unique<RzGridVolumeLookup>(gctx, *world, config, *logger));
    BOOST_CHECK(output.str().empty());
    const auto before = calls;
    BOOST_CHECK(world->lowestTrackingVolume(gctx, {15., 0., 0.}) == expected);
    BOOST_CHECK_EQUAL(calls, before + 1);
    BOOST_CHECK(world->lowestTrackingVolume(gctx, {40., 0., 0.}) ==
                &world->volumes().at(0));
    BOOST_CHECK_EQUAL(calls, before + 2);
    BOOST_CHECK(world->lowestTrackingVolume(gctx, {70., 0., 0.}) ==
                world.get());
    BOOST_CHECK_EQUAL(calls, before + 2);
  }
}

BOOST_AUTO_TEST_CASE(StaticNodeCanAccelerateWithoutWorldPolicy) {
  Blueprint blueprint{Blueprint::Config{}};
  auto region =
      std::make_shared<StaticBlueprintNode>(std::make_unique<TrackingVolume>(
          Transform3::Identity(),
          std::make_shared<CylinderVolumeBounds>(0., 20., 20.), "region"));
  region->addChild(
      std::make_shared<StaticBlueprintNode>(std::make_unique<TrackingVolume>(
          Transform3::Identity(),
          std::make_shared<CylinderVolumeBounds>(0., 10., 10.), "leaf")));
  region->setVolumeLookupFactory([](const GeometryContext& context,
                                    const TrackingVolume& volume,
                                    const Logger& logger) {
    return std::make_unique<RzGridVolumeLookup>(
        context, volume, RzGridVolumeLookup::Config{}, logger);
  });
  blueprint.addChild(region);
  auto geometry = blueprint.construct(BlueprintOptions{}, gctx);
  BOOST_CHECK(geometry->highestTrackingVolume()->volumeLookup() == nullptr);
  const auto& builtRegion = geometry->highestTrackingVolume()->volumes().at(0);
  BOOST_CHECK(dynamic_cast<const RzGridVolumeLookup*>(
                  builtRegion.volumeLookup()) != nullptr);
  BOOST_CHECK(builtRegion.volumes().at(0).volumeLookup() == nullptr);
  const auto* found =
      *geometry->resolveLowestTrackingVolume(gctx, {5., 0., 0.});
  BOOST_REQUIRE(found != nullptr);
  BOOST_CHECK_EQUAL(found->volumeName(), "leaf");
}

BOOST_AUTO_TEST_CASE(LinearDiagnosticsRespectSubtreeBoundaries) {
  // These assertions deliberately exercise warning diagnostics.
#ifdef ACTS_ENABLE_LOG_FAILURE_THRESHOLD
  Logging::ScopedFailureThreshold threshold{Logging::ERROR};
#endif
  auto world = std::make_unique<TrackingVolume>(
      Transform3::Identity(),
      std::make_shared<CylinderVolumeBounds>(0., 100., 100.));
  auto region = std::make_unique<TrackingVolume>(
      Transform3::Identity(),
      std::make_shared<CylinderVolumeBounds>(0., 50., 50.));
  for (int i = 0; i < 65; ++i) {
    region->addVolume(std::make_unique<TrackingVolume>(
        Transform3::Identity(),
        std::make_shared<CylinderVolumeBounds>(0., 10., 10.)));
  }
  auto* regionPtr = region.get();
  world->addVolume(std::move(region));
  std::ostringstream output;
  auto logger = getDefaultLogger("LinearScopeTest", Logging::WARNING, &output);
  LinearVolumeLookup scanned(*world, LinearVolumeLookup::Config{}, *logger);
  BOOST_CHECK(output.str().find("scans up to 65 siblings") !=
              std::string::npos);
  output.str("");
  regionPtr->setVolumeLookup(std::make_unique<RzGridVolumeLookup>(
      gctx, *regionPtr, RzGridVolumeLookup::Config{}));
  LinearVolumeLookup delegated(*world, LinearVolumeLookup::Config{}, *logger);
  BOOST_CHECK(output.str().empty());
}

BOOST_AUTO_TEST_CASE(RzLookupCoincidentBounds) {
  auto world = std::make_shared<TrackingVolume>(
      Transform3::Identity(),
      std::make_shared<CylinderVolumeBounds>(0., 100., 100.));
  // Coincident rectangles share a cell and must retain sibling priority.
  for (int i = 0; i < 200; ++i) {
    world->addVolume(std::make_unique<TrackingVolume>(
        Transform3::Identity(),
        std::make_shared<CylinderVolumeBounds>(0., 10., 10.)));
  }
  {
    RzGridVolumeLookup lookup(gctx, *world, RzGridVolumeLookup::Config{});
    BOOST_CHECK(lookup.lowestTrackingVolume(gctx, {5., 0., 0.}, 0.) ==
                &world->volumes().at(0));
    BOOST_CHECK(lookup.lowestTrackingVolume(gctx, {15., 0., 0.}, 0.) ==
                world.get());
  }
}

BOOST_AUTO_TEST_CASE(RzLookupsRejectGen1Storage) {
  auto makeGen1 = [] {
    auto dense = std::make_shared<TrackingVolume>(
        Transform3::Identity(),
        std::make_shared<CylinderVolumeBounds>(0., 5., 5.));
    return std::make_unique<TrackingVolume>(
        Transform3::Identity(),
        std::make_shared<CylinderVolumeBounds>(0., 20., 20.), nullptr, nullptr,
        nullptr, MutableTrackingVolumeVector{dense}, "gen1");
  };
  {
    auto gen1 = makeGen1();
    BOOST_CHECK_THROW(
        (RzGridVolumeLookup(gctx, *gen1, RzGridVolumeLookup::Config{})),
        std::invalid_argument);
    // Implicit and explicit linear lookup still retain Gen1 semantics.
    const Vector3 point(2., 0., 0.);
    const auto* dense = gen1->denseVolumes().front().get();
    BOOST_CHECK(gen1->lowestTrackingVolume(gctx, point) == dense);
    gen1->setVolumeLookup(std::make_unique<LinearVolumeLookup>(
        *gen1, LinearVolumeLookup::Config{}));
    BOOST_CHECK(gen1->lowestTrackingVolume(gctx, point) == dense);

    auto world = makeWorld();
    auto* gen1Ptr = gen1.get();
    world->addVolume(std::move(gen1));
    // Also reject Gen1 storage at a descendant, even at a lookup boundary.
    BOOST_CHECK_THROW(
        (RzGridVolumeLookup(gctx, *world, RzGridVolumeLookup::Config{})),
        std::invalid_argument);
    gen1Ptr->setVolumeLookup(nullptr);
    BOOST_CHECK_THROW(
        (RzGridVolumeLookup(gctx, *world, RzGridVolumeLookup::Config{})),
        std::invalid_argument);
  }
}

BOOST_AUTO_TEST_CASE(RzGridNearIdentityRotation) {
  auto world = std::make_shared<TrackingVolume>(
      Transform3::Identity(),
      std::make_shared<CylinderVolumeBounds>(0., 100., 2000.));
  Transform3 transform = Transform3::Identity();
  transform.linear() = AngleAxis3(5e-13, Vector3::UnitY()).toRotationMatrix();
  auto child = std::make_unique<TrackingVolume>(
      transform, std::make_shared<CylinderVolumeBounds>(5., 20., 1000.));
  const auto* expected = child.get();
  world->addVolume(std::move(child));
  std::ostringstream output;
  auto logger = getDefaultLogger("RotationTest", Logging::WARNING, &output);
  RzGridVolumeLookup::Config config;
  config.fallbackWarningFraction = 0.;
  RzGridVolumeLookup lookup(gctx, *world, config, *logger);
  BOOST_CHECK(output.str().empty());
  for (double z : {-1000., 0., 1000.}) {
    for (double r : {5., 20.}) {
      const Vector3 point = transform * Vector3(r, 0., z);
      BOOST_CHECK(lookup.lowestTrackingVolume(gctx, point, 1e-6) == expected);
    }
  }
}

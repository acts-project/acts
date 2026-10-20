// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/Units.hpp"
#include "Acts/Geometry/CylinderPortalShell.hpp"
#include "Acts/Geometry/CylinderVolumeBounds.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/TrackingVolume.hpp"
#include "Acts/Navigation/CylinderNavigationPolicy.hpp"
#include "Acts/Navigation/NavigationStream.hpp"
#include "Acts/Navigation/TryAllNavigationPolicy.hpp"
#include "Acts/Surfaces/CylinderSurface.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "ActsTests/CommonHelpers/BenchmarkTools.hpp"

#include <array>
#include <chrono>
#include <cmath>
#include <iostream>
#include <memory>
#include <numbers>
#include <random>
#include <string_view>
#include <vector>

using namespace Acts;
using namespace Acts::UnitLiterals;

namespace ActsTests {
namespace {

constexpr double rMin = 100_mm;
constexpr double rMax = 400_mm;
constexpr double halfLengthZ = 300_mm;
constexpr std::size_t nInputs = 4096;
constexpr std::size_t nRuns = 2000;
constexpr auto warmupTime = std::chrono::milliseconds(500);

struct Scenario {
  std::string_view name;
  Transform3 transform;
};

std::array<Scenario, 3> scenarios() {
  Transform3 translated = Transform3::Identity();
  translated.translation() = Vector3(30_mm, -20_mm, 50_mm);
  Transform3 rotated = translated;
  rotated.rotate(AngleAxis3(0.35, Vector3(1., 2., 3.).normalized()));
  return {{{"aligned", Transform3::Identity()},
           {"translated", translated},
           {"rotated", rotated}}};
}

std::vector<NavigationArguments> makeInputs(const Transform3& transform) {
  std::mt19937 generator(42);
  std::uniform_real_distribution<double> uniform(0., 1.);
  std::vector<NavigationArguments> inputs;
  inputs.reserve(nInputs);
  for (std::size_t i = 0; i < nInputs; ++i) {
    // Sample points uniformly in the shell volume. The same local rays are
    // used for all transforms, keeping the geometric workload comparable.
    const double radius = std::sqrt(
        rMin * rMin + uniform(generator) * (rMax * rMax - rMin * rMin));
    const double phi = 2. * std::numbers::pi * uniform(generator);
    const double z = halfLengthZ * (2. * uniform(generator) - 1.);
    const Vector3 position(radius * std::cos(phi), radius * std::sin(phi), z);

    // Mix outward tracks from the interaction region with isotropic tracks
    // to cover inner/outer cylinder hits, misses, and exits through discs.
    const double dirPhi = 2. * std::numbers::pi * uniform(generator);
    const double cosTheta = 2. * uniform(generator) - 1.;
    const double sinTheta = std::sqrt(1. - cosTheta * cosTheta);
    Vector3 direction(sinTheta * std::cos(dirPhi), sinTheta * std::sin(dirPhi),
                      cosTheta);
    if (i % 4 == 0) {
      direction = position.normalized();
    }
    inputs.push_back({.position = transform * position,
                      .direction = transform.linear() * direction});
  }
  return inputs;
}

template <bool InitializeStream, typename Policy>
MicroBenchmarkResult benchmarkPolicy(
    const GeometryContext& gctx, const Policy& policy,
    const std::vector<NavigationArguments>& inputs, const Logger& logger) {
  NavigationStream stream;
  stream.reserve(8);
  AppendOnlyNavigationStream appendStream(stream);
  NavigationPolicyStateManager stateManager;
  policy.createState(gctx, inputs.front(), stateManager, logger);
  auto state = stateManager.currentState();
  auto iteration = [&](const NavigationArguments& args) {
    stream.reset();
    policy.initializeCandidates(gctx, args, state, appendStream, logger);
    if constexpr (InitializeStream) {
      // Both policies add each portal only once, just as when the navigator
      // initializes candidates for a new volume.
      stream.initialize(gctx, {args.position, args.direction}, logger,
                        s_onSurfaceTolerance, true);
    }
    return stream.candidates().size();
  };

  // Prewarm every query so both the main and additional-intersection buffers
  // have reached their maximum capacity before any timed measurements.
  for (const auto& input : inputs) {
    iteration(input);
  }
  return microBenchmark(iteration, inputs, nRuns, warmupTime);
}

}  // namespace

BOOST_AUTO_TEST_CASE(cylinder_intersections) {
  const auto gctx = GeometryContext::dangerouslyDefaultConstruct();
  for (const auto& [name, transform] : scenarios()) {
    const auto inputs = makeInputs(transform);
    const auto inner =
        Surface::makeShared<CylinderSurface>(transform, rMin, halfLengthZ);
    const auto outer =
        Surface::makeShared<CylinderSurface>(transform, rMax, halfLengthZ);
    for (bool bounded : {true, false}) {
      const auto tolerance =
          bounded ? BoundaryTolerance::None() : BoundaryTolerance::Infinite();
      const auto result = microBenchmark(
          [&](const NavigationArguments& args) {
            return std::array{inner->intersect(gctx, args.position,
                                               args.direction, tolerance),
                              outer->intersect(gctx, args.position,
                                               args.direction, tolerance)};
          },
          inputs, nRuns, warmupTime);
      std::cout << name << "/" << (bounded ? "None" : "Infinite")
                << " (inner+outer intersection pair): " << result << '\n';
    }
  }
}

BOOST_AUTO_TEST_CASE(cylinder_navigation) {
  const auto gctx = GeometryContext::dangerouslyDefaultConstruct();
  const auto logger =
      getDefaultLogger("CylinderNavigationBenchmark", Logging::WARNING);
  for (const auto& [name, transform] : scenarios()) {
    const auto inputs = makeInputs(transform);
    const auto bounds =
        std::make_shared<CylinderVolumeBounds>(rMin, rMax, halfLengthZ);
    TrackingVolume volume(transform, bounds, "Cylinder");
    SingleCylinderPortalShell shell(gctx, volume);
    shell.applyToVolume();
    const CylinderNavigationPolicy cylinder(gctx, volume, *logger);
    const TryAllNavigationPolicy tryAll(
        gctx, volume, *logger,
        {.portals = true, .sensitives = false, .passives = false});
    auto benchmark = [&](const auto& policy, std::string_view policyName) {
      std::cout << name << "/" << policyName << " (candidates): "
                << benchmarkPolicy<false>(gctx, policy, inputs, *logger)
                << '\n';
      std::cout << name << "/" << policyName << " (candidates+stream): "
                << benchmarkPolicy<true>(gctx, policy, inputs, *logger) << '\n';
    };
    benchmark(cylinder, "Cylinder");
    benchmark(tryAll, "TryAll");
  }
}

}  // namespace ActsTests

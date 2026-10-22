// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Geometry/CylinderVolumeBounds.hpp"
#include "Acts/Geometry/RzGridVolumeLookup.hpp"
#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/Geometry/TrackingVolume.hpp"
#include "ActsPlugins/Json/TrackingGeometryJsonConverter.hpp"

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <limits>
#include <memory>
#include <numbers>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
using Clock = std::chrono::steady_clock;
using Volume = Acts::TrackingVolume;
using Vec = Acts::Vector3;
constexpr double tolerance = 1e-6;
template <typename F>
double elapsed(F&& f) {
  const auto start = Clock::now();
  f();
  return std::chrono::duration<double>(Clock::now() - start).count();
}

// Output checksum makes all timed results observable; validation is untimed.
template <typename F>
void measure(const std::string& workload, const std::string& method,
             const std::vector<Vec>& points, std::size_t repeats, F&& lookup) {
  std::vector<double> times;
  std::uint64_t checksum = 0;
  for (std::size_t r = 0; r < repeats; ++r) {
    const double seconds = elapsed([&] {
      for (const auto& p : points) {
        const auto* v = lookup(p);
        checksum += v == nullptr ? 0 : v->geometryId().value();
      }
    });
    times.push_back(seconds * 1e9 / points.size());
  }
  std::sort(times.begin(), times.end());
  std::cout << workload << ',' << method << ',' << points.size() << ','
            << times[times.size() / 2] << ',' << times.front() << ','
            << checksum << '\n';
}
}  // namespace

int main(int argc, char** argv) try {
  if (argc < 2 || argc > 4) {
    std::cerr << "Usage: " << argv[0]
              << " geometry.json [points=100000] [repeats=7]\n";
    return 1;
  }
  const auto count = argc > 2 ? std::stoull(argv[2]) : 100000;
  const auto repeats = argc > 3 ? std::stoull(argv[3]) : 7;
  if (count == 0 || repeats == 0) {
    throw std::invalid_argument("points and repeats must be positive");
  }
  const auto gctx = Acts::GeometryContext::dangerouslyDefaultConstruct();
  Acts::TrackingGeometryJsonConverter converter(
      Acts::TrackingGeometryJsonConverter::Config::defaultConfig(),
      Acts::getDefaultLogger("VolumeLookupBenchmark", Acts::Logging::WARNING,
                             &std::cerr));
  std::ifstream input(argv[1]);
  if (!input) {
    throw std::runtime_error("Cannot open geometry file");
  }
  auto payload = nlohmann::json::parse(input);
  // Older exports include two bevel angles, removed from CylinderVolumeBounds.
  // Only zero bevels have an equivalent representation in the current schema.
  std::size_t normalized = 0;
  for (auto& volume : payload.at("volumes")) {
    auto& bounds = volume.at("bounds");
    auto& values = bounds.at("values");
    if (bounds.at("kind") == "Cylinder" && values.size() == 7) {
      if (values.at(5).get<double>() != 0 || values.at(6).get<double>() != 0) {
        throw std::runtime_error(
            "Cannot normalize nonzero legacy cylinder bevels");
      }
      values.erase(values.begin() + 5, values.end());
      ++normalized;
    }
  }
  std::cerr << "Normalized zero-bevel legacy cylinders: " << normalized << '\n';
  auto geometry = converter.fromJson(gctx, payload);
  const auto& world = *geometry->highestTrackingVolume();
  const auto baseline = [&](const Vec& p) {
    return world.lowestTrackingVolume(gctx, p, tolerance);
  };
  auto productionGrid = converter.fromJson(gctx, payload);
  const Acts::RzGridVolumeLookup::Config config;
  std::cerr << "Grid build ms: "
            << 1000 * elapsed([&] {
                 productionGrid->highestTrackingVolume()->setVolumeLookup(
                     std::make_unique<Acts::RzGridVolumeLookup>(
                         gctx, *productionGrid->highestTrackingVolume(),
                         config));
               })
            << '\n';
  const auto sameVolume = [](const Volume* lhs, const Volume* rhs) {
    return lhs == nullptr || rhs == nullptr
               ? lhs == rhs
               : lhs->geometryId() == rhs->geometryId();
  };
  const auto productionGridLookup = [&](const Vec& p) {
    return *productionGrid->resolveLowestTrackingVolume(gctx, p, std::nullopt,
                                                        nullptr, tolerance);
  };
  // Exercise partition edges themselves, including adjacent representable
  // coordinates. These checks are untimed and independent of random workloads.
  std::size_t edgeChecks = 0;
  auto checkEdge = [&](double r, double z) {
    if (r < 0) {
      return;
    }
    for (double direction : {-std::numeric_limits<double>::infinity(), 0.,
                             std::numeric_limits<double>::infinity()}) {
      Vec p(r, 0, z);
      if (direction != 0) {
        p.x() = std::nextafter(p.x(), direction);
        p.z() = std::nextafter(p.z(), direction);
      }
      const auto* expected = baseline(p);
      if (!sameVolume(productionGridLookup(p), expected)) {
        throw std::runtime_error("Lookup mismatch at partition edge");
      }
      ++edgeChecks;
    }
  };
  // Reproduce the grid envelope edges of this unrotated-cylinder export.
  std::array<std::vector<double>, 2> edges;
  world.visitVolumes([&](const Volume* volume) {
    if (volume == &world) {
      return;
    }
    const auto& bounds =
        dynamic_cast<const Acts::CylinderVolumeBounds&>(volume->volumeBounds());
    const auto& transform = volume->localToGlobalTransform(gctx);
    if (!transform.linear().isIdentity(0)) {
      throw std::runtime_error("Edge validation requires unrotated cylinders");
    }
    const double margin = config.maxTolerance + 1e-9;
    const double offset = transform.translation().head<2>().norm();
    edges[0].push_back(bounds.get(Acts::CylinderVolumeBounds::eMinR) - offset -
                       margin);
    edges[0].push_back(bounds.get(Acts::CylinderVolumeBounds::eMaxR) + offset +
                       margin);
    edges[1].push_back(transform.translation().z() -
                       bounds.get(Acts::CylinderVolumeBounds::eHalfLengthZ) -
                       margin);
    edges[1].push_back(transform.translation().z() +
                       bounds.get(Acts::CylinderVolumeBounds::eHalfLengthZ) +
                       margin);
  });
  for (auto& axis : edges) {
    std::sort(axis.begin(), axis.end());
    axis.erase(std::unique(axis.begin(), axis.end()), axis.end());
  }
  for (double r : edges[0]) {
    for (double z : edges[1]) {
      checkEdge(r, z);
    }
  }
  std::cerr << "Validated partition-edge queries: " << edgeChecks << '\n';
  std::vector<const Volume*> volumes;
  world.visitVolumes([&](const Volume* v) { volumes.push_back(v); });
  std::cerr << "Volumes: " << volumes.size()
            << "; seed: 42; tolerance: " << tolerance << '\n';
  std::mt19937_64 random(42);
  std::uniform_real_distribution<double> uniform(0, 1);
  const auto bounds = world.boundingBox(gctx);
  std::cout << "workload,method,points,median_ns,min_ns,checksum\n";
  for (const std::string workload :
       {"world-box", "volume-balanced", "boundaries", "outside"}) {
    std::vector<Vec> points;
    points.reserve(count);
    while (points.size() < count) {
      Vec p;
      if (workload == "world-box" || workload == "outside") {
        for (int axis = 0; axis < 3; ++axis) {
          p[axis] = bounds.min()[axis] +
                    uniform(random) * (bounds.max()[axis] - bounds.min()[axis]);
        }
        if (workload == "outside") {
          p.z() = bounds.max().z() + 1 + uniform(random) * 100;
        }
      } else {
        const auto& v = *volumes[points.size() % volumes.size()];
        const auto& b =
            dynamic_cast<const Acts::CylinderVolumeBounds&>(v.volumeBounds());
        const double rmin = b.get(Acts::CylinderVolumeBounds::eMinR);
        const double rmax = b.get(Acts::CylinderVolumeBounds::eMaxR);
        const double hz = b.get(Acts::CylinderVolumeBounds::eHalfLengthZ);
        double r = std::sqrt(rmin * rmin +
                             uniform(random) * (rmax * rmax - rmin * rmin));
        double z = (2 * uniform(random) - 1) * hz;
        if (workload == "boundaries") {
          const auto sample = points.size() / volumes.size();
          const double delta =
              std::array{-2., -0.5, 0., 0.5, 2.}[(sample / 4) % 5] * tolerance;
          switch (sample % 4) {
            case 0:
              r = std::max(0., rmin + delta);
              break;
            case 1:
              r = rmax + delta;
              break;
            case 2:
              z = -hz + delta;
              break;
            case 3:
              z = hz + delta;
              break;
          }
        }
        const double phi = 2 * std::numbers::pi * uniform(random);
        p = v.localToGlobalTransform(gctx) *
            Vec(r * std::cos(phi), r * std::sin(phi), z);
      }
      points.push_back(p);
    }
    std::shuffle(points.begin(), points.end(), random);
    std::size_t misses = 0;
    for (const auto& p : points) {
      const auto* expected = baseline(p);
      if (expected == nullptr) {
        ++misses;
      }
      if (!sameVolume(productionGridLookup(p), expected)) {
        throw std::runtime_error("Lookup mismatch in " + workload);
      }
    }
    std::cerr << workload << ": validated " << points.size()
              << ", outside world " << misses << '\n';
    measure(workload, "linear", points, repeats, baseline);
    measure(workload, "production-rz-grid", points, repeats,
            productionGridLookup);
  }
  return 0;
} catch (const std::exception& e) {
  std::cerr << e.what() << '\n';
  return 1;
}

// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Material/IntersectionMaterialAssigner.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/Surfaces/RectangleBounds.hpp"

#include <chrono>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <random>
#include <sstream>

#include "MaterialAssignment/IndexedMaterialAssigner.hpp"

using namespace Acts;
using ActsTests::IndexedMaterialAssigner;
using Clock = std::chrono::steady_clock;

namespace {
struct RayInput {
  Vector3 position, direction;
};

std::vector<std::shared_ptr<Surface>> readSurfaces(const std::string& path) {
  std::ifstream input(path);
  if (!input) {
    throw std::runtime_error("Cannot open surface fixture: " + path);
  }
  std::vector<std::shared_ptr<Surface>> surfaces;
  std::string line;
  while (std::getline(input, line)) {
    if (line.empty() || line.front() == '#') {
      continue;
    }
    std::istringstream row(line);
    std::string type;
    std::uint64_t id;
    std::array<double, 4> bounds{};
    Transform3 transform = Transform3::Identity();
    row >> type >> id;
    for (auto& b : bounds) {
      row >> b;
    }
    for (int i = 0; i < 3; ++i) {
      for (int j = 0; j < 3; ++j) {
        row >> transform.linear()(i, j);
      }
    }
    for (int i = 0; i < 3; ++i) {
      row >> transform.translation()[i];
    }
    if (!row) {
      throw std::runtime_error("Malformed surface row");
    }
    std::shared_ptr<Surface> surface;
    if (type == "cylinder") {
      surface = Surface::makeShared<CylinderSurface>(
          transform, std::make_shared<CylinderBounds>(bounds));
    } else if (type == "disc") {
      surface = Surface::makeShared<DiscSurface>(
          transform, std::make_shared<RadialBounds>(bounds));
    } else {
      throw std::runtime_error("Unknown surface type: " + type);
    }
    surface->assignGeometryId(GeometryIdentifier(id));
    surfaces.push_back(std::move(surface));
  }
  if (surfaces.empty()) {
    throw std::runtime_error("Empty surface fixture");
  }
  return surfaces;
}

void compare(const IAssignmentFinder& baseline, const IAssignmentFinder& test,
             const GeometryContext& gctx, const std::vector<RayInput>& rays) {
  const MagneticFieldContext mctx;
  for (std::size_t i = 0; i < rays.size(); ++i) {
    const auto& ray = rays[i];
    auto a =
        baseline.assignmentCandidates(gctx, mctx, ray.position, ray.direction);
    auto b = test.assignmentCandidates(gctx, mctx, ray.position, ray.direction);
    bool equal =
        a.first.size() == b.first.size() && a.second.size() == b.second.size();
    if (equal) {
      for (std::size_t j = 0; j < a.first.size(); ++j) {
        equal &= a.first[j].surface == b.first[j].surface &&
                 (a.first[j].position - b.first[j].position).norm() < 1e-9 &&
                 a.first[j].direction == b.first[j].direction;
      }
    }
    if (!equal) {
      std::ostringstream message;
      message << "Mismatch at ray " << i
              << " origin=" << ray.position.transpose()
              << " direction=" << ray.direction.transpose()
              << " hits=" << a.first.size() << "/" << b.first.size();
      throw std::runtime_error(message.str());
    }
  }
}

std::size_t checksum = 0;
double measure(const IAssignmentFinder& finder, const GeometryContext& gctx,
               const std::vector<RayInput>& rays) {
  const MagneticFieldContext mctx;
  const auto start = Clock::now();
  std::size_t sum = 0;
  for (const auto& ray : rays) {
    auto hits =
        finder.assignmentCandidates(gctx, mctx, ray.position, ray.direction);
    sum += hits.first.size();
  }
  checksum += sum;
  return std::chrono::duration<double, std::micro>(Clock::now() - start)
             .count() /
         static_cast<double>(rays.size());
}
// Compare the complete two-intersection API, not only closestForward().
void cylinderExperiment(const GeometryContext& gctx, std::size_t count,
                        std::size_t repeats) {
  auto cylinder = Surface::makeShared<CylinderSurface>(
      Transform3::Identity(), std::make_shared<CylinderBounds>(100, 400));
  std::mt19937_64 rng(314159);
  std::uniform_real_distribution<double> u(-1., 1.);
  std::vector<RayInput> onAxis, offAxis;
  for (std::size_t i = 0; i < count; ++i) {
    Vector3 d(u(rng), u(rng), u(rng));
    d.normalize();
    onAxis.push_back({Vector3(0, 0, 100 * u(rng)), d});
    offAxis.push_back({Vector3(10 * u(rng), 10 * u(rng), 100 * u(rng)), d});
  }
  auto validate = [&](const CylinderSurface& surface, const RayInput& ray,
                      const BoundaryTolerance& boundary, double tolerance) {
    auto a = surface.intersect(gctx, ray.position, ray.direction, boundary,
                               tolerance);
    auto b = ActsTests::intersectOnAxisCylinder(
        surface, gctx, ray.position, ray.direction, boundary, tolerance);
    if (a.size() != b.size()) {
      throw std::runtime_error("Cylinder root count mismatch");
    }
    for (unsigned int i = 0; i < a.size(); ++i) {
      auto equalNumber = [](double x, double y) {
        return (std::isnan(x) && std::isnan(y)) || x == y ||
               (std::isfinite(x) && std::isfinite(y) &&
                std::abs(x - y) <=
                    1e-10 * std::max({1., std::abs(x), std::abs(y)}));
      };
      bool equal = a[i].status() == b[i].status() &&
                   equalNumber(a[i].pathLength(), b[i].pathLength());
      for (int axis = 0; axis < 3; ++axis) {
        equal &= equalNumber(a[i].position()[axis], b[i].position()[axis]);
      }
      if (!equal) {
        throw std::runtime_error("Cylinder root/status mismatch");
      }
    }
  };
  for (const auto& boundary :
       {BoundaryTolerance::None(), BoundaryTolerance::Infinite(),
        BoundaryTolerance::AbsoluteEuclidean(0.01)}) {
    for (const auto& ray : onAxis) {
      validate(*cylinder, ray, boundary, s_onSurfaceTolerance);
    }
    for (const auto& ray : offAxis) {
      validate(*cylinder, ray, boundary, s_onSurfaceTolerance);
    }
    for (double z : {0., 400., 400. + s_onSurfaceTolerance, 401.}) {
      for (const Vector3& direction :
           {Vector3(1, 0, 0), Vector3(0, 0, 1), Vector3(0, 0, -1),
            Vector3(2, 0, 0), Vector3(1e-14, 0, 1)}) {
        for (double tolerance : {s_onSurfaceTolerance, 200.}) {
          validate(*cylinder, {Vector3(0, 0, z), direction}, boundary,
                   tolerance);
        }
      }
    }
  }
  // A translated cylinder exercises the same path;
  // partial-phi cylinders exercise the fallback.
  auto transformed = Surface::makeShared<CylinderSurface>(
      Transform3(Translation3(20, 30, 40)),
      std::make_shared<CylinderBounds>(100, 400));
  auto sector = Surface::makeShared<CylinderSurface>(
      Transform3::Identity(), std::make_shared<CylinderBounds>(100, 400, 0.5));
  Transform3 rotatedPlacement = Transform3::Identity();
  rotatedPlacement.linear() << 0, 0, 1, 1, 0, 0, 0, 1, 0;
  auto rotated = Surface::makeShared<CylinderSurface>(
      rotatedPlacement, std::make_shared<CylinderBounds>(100, 400));
  for (const auto& ray : onAxis) {
    validate(*rotated,
             {rotatedPlacement * ray.position,
              rotatedPlacement.linear() * ray.direction},
             BoundaryTolerance::None(), s_onSurfaceTolerance);
    validate(*transformed, {ray.position + Vector3(20, 30, 40), ray.direction},
             BoundaryTolerance::None(), s_onSurfaceTolerance);
    validate(*sector, ray, BoundaryTolerance::None(), s_onSurfaceTolerance);
  }
  std::cout << "cylinder_two_root_checks=passed\n";
  double sum = 0;
  for (bool displaced : {false, true}) {
    const auto& rays = displaced ? offAxis : onAxis;
    std::array<std::vector<double>, 2> timings;
    for (std::size_t repeat = 0; repeat <= repeats; ++repeat) {
      for (std::size_t j = 0; j < 2; ++j) {
        auto k = (j + repeat) % 2;
        const auto start = Clock::now();
        for (const auto& ray : rays) {
          auto hits =
              k == 0 ? cylinder->intersect(gctx, ray.position, ray.direction,
                                           BoundaryTolerance::None())
                     : ActsTests::intersectOnAxisCylinder(
                           *cylinder, gctx, ray.position, ray.direction,
                           BoundaryTolerance::None());
          sum += hits[1].pathLength();
          sum += static_cast<int>(hits[0].status()) +
                 static_cast<int>(hits[1].status());
        }
        const double ns =
            std::chrono::duration<double, std::nano>(Clock::now() - start)
                .count() /
            static_cast<double>(rays.size());
        if (repeat != 0) {
          timings[k].push_back(ns);
        }
      }
    }
    for (std::size_t k = 0; k < 2; ++k) {
      std::ranges::sort(timings[k]);
      std::cout << "cylinder " << (displaced ? "off-axis" : "on-axis") << ' '
                << (k == 0 ? "baseline" : "fast")
                << " ns/intersect=" << timings[k][repeats / 2]
                << " min=" << timings[k].front() << " max=" << timings[k].back()
                << '\n';
    }
  }
  std::cout << "cylinder_checksum=" << sum << '\n';
}

}  // namespace

int main(int argc, char** argv) {
  try {
    if (argc < 2 || argc > 4) {
      std::cerr << "Usage: " << argv[0]
                << " surfaces.txt [rays=100000] [repeats=7]\n";
      return 1;
    }
    const std::size_t count = argc > 2 ? std::stoul(argv[2]) : 100000;
    const std::size_t repeats = argc > 3 ? std::stoul(argv[3]) : 7;
    if (count == 0 || repeats == 0) {
      throw std::runtime_error("rays and repeats must be positive");
    }
    auto owned = readSurfaces(argv[1]);
    std::vector<const Surface*> surfaces;
    for (const auto& s : owned) {
      surfaces.push_back(s.get());
    }
    const auto gctx = GeometryContext::dangerouslyDefaultConstruct();
    IntersectionMaterialAssigner baseline(
        {surfaces, {}}, getDefaultLogger("baseline", Logging::ERROR));
    std::cout << "surfaces=" << surfaces.size() << " rays=" << count
              << " repeats=" << repeats << " seed=20261015\n";
    std::vector<std::unique_ptr<IndexedMaterialAssigner>> alternatives;
    std::vector<const char*> names{"baseline",
                                   "octree",
                                   "rz-tree",
                                   "rz-linear",
                                   "angular",
                                   "brute-control",
                                   "brute-fast-cylinder",
                                   "angular-fast-cylinder"};
    for (auto mode : {IndexedMaterialAssigner::Mode::Octree,
                      IndexedMaterialAssigner::Mode::RZ,
                      IndexedMaterialAssigner::Mode::LinearRZ,
                      IndexedMaterialAssigner::Mode::Angular,
                      IndexedMaterialAssigner::Mode::Brute,
                      IndexedMaterialAssigner::Mode::Brute,
                      IndexedMaterialAssigner::Mode::Angular}) {
      auto start = Clock::now();
      alternatives.push_back(std::make_unique<IndexedMaterialAssigner>(
          gctx, surfaces, mode, alternatives.size() >= 5));
      std::cout << names[alternatives.size()] << " build_us="
                << std::chrono::duration<double, std::micro>(Clock::now() -
                                                             start)
                       .count()
                << " fallback=" << alternatives.back()->fallbackSize() << '\n';
    }
    std::mt19937_64 rng(20261015);
    std::uniform_real_distribution<double> uniform(0., 1.);
    std::vector<RayInput> edgeRays;
    for (const auto& s : owned) {
      const auto& trf = s->localToGlobalTransform(gctx);
      const double z = trf.translation().z();
      double radius = s->bounds().values()[0];
      if (auto* cb = dynamic_cast<const CylinderBounds*>(&s->bounds())) {
        for (double sign : {-1., 1.}) {
          for (double eps : {-1e-3, 0., 1e-3}) {
            edgeRays.push_back(
                {Vector3(
                     0, 0,
                     z + sign * cb->get(CylinderBounds::eHalfLengthZ) + eps),
                 Vector3::UnitX()});
          }
        }
      }
      for (const Vector3& p :
           {Vector3(0, 0, z), Vector3(radius, 0, z), Vector3(radius, -100, z),
            Vector3(radius, 0, z - 100)}) {
        for (const Vector3& d :
             {Vector3::UnitX().eval(), Vector3::UnitY().eval(),
              Vector3::UnitZ().eval(), (-Vector3::UnitZ()).eval(),
              Vector3(1, 0, 1e-14).normalized()}) {
          edgeRays.push_back({p, d});
        }
      }
    }
    // Exercise angular endpoints and their immediate floating-point neighbours.
    for (const auto& surface : owned) {
      const double z = surface->localToGlobalTransform(gctx).translation().z();
      const auto values = surface->bounds().values();
      const bool cylinder = surface->type() == Surface::Cylinder;
      const double pad = 2 * s_onSurfaceTolerance;
      const double rmin = std::max(0., values[0] - pad);
      const double rmax = (cylinder ? values[0] : values[1]) + pad;
      const double zmin = z - (cylinder ? values[1] : 0.) - pad;
      const double zmax = z + (cylinder ? values[1] : 0.) + pad;
      if (rmin == 0) {
        continue;
      }
      for (double q : {zmin / rmin, zmin / rmax, zmax / rmin, zmax / rmax}) {
        for (double slope :
             {std::nextafter(q, -std::numeric_limits<double>::infinity()), q,
              std::nextafter(q, std::numeric_limits<double>::infinity())}) {
          edgeRays.push_back({Vector3::Zero(), Vector3(1, 0, slope)});
        }
      }
    }
    for (const auto& alternative : alternatives) {
      compare(baseline, *alternative, gctx, edgeRays);
    }
    std::cout << "edge_checks=" << edgeRays.size() << " passed\n";
    for (const std::string sample : {"ip", "beamspot", "displaced"}) {
      std::vector<RayInput> rays;
      rays.reserve(count);
      for (std::size_t i = 0; i < count; ++i) {
        const double eta = 8 * uniform(rng) - 4;
        const double phi = 2 * std::numbers::pi * uniform(rng);
        Vector3 p = Vector3::Zero();
        if (sample == "beamspot") {
          p = Vector3(0.1 * (uniform(rng) - 0.5), 0.1 * (uniform(rng) - 0.5),
                      100 * (uniform(rng) - 0.5));
        } else if (sample == "displaced") {
          p = Vector3(2400 * (uniform(rng) - 0.5), 2400 * (uniform(rng) - 0.5),
                      8000 * (uniform(rng) - 0.5));
        }
        rays.push_back(
            {p, Vector3(std::cos(phi) / std::cosh(eta),
                        std::sin(phi) / std::cosh(eta), std::tanh(eta))});
      }
      for (const auto& alternative : alternatives) {
        compare(baseline, *alternative, gctx, rays);
      }
      std::vector<const IAssignmentFinder*> finders{&baseline};
      for (const auto& alternative : alternatives) {
        finders.push_back(alternative.get());
      }
      std::vector<std::vector<double>> timings(finders.size());
      for (auto* finder : finders) {
        measure(*finder, gctx, rays);
      }
      // Rotate order each repetition to reduce systematic warmup/order bias.
      for (std::size_t repeat = 0; repeat < repeats; ++repeat) {
        for (std::size_t j = 0; j < finders.size(); ++j) {
          const auto k = (j + repeat) % finders.size();
          timings[k].push_back(measure(*finders[k], gctx, rays));
        }
      }
      std::cout << sample << " equivalence=passed\n";
      for (std::size_t k = 0; k < timings.size(); ++k) {
        std::ranges::sort(timings[k]);
        double candidates = static_cast<double>(surfaces.size());
        if (k > 0) {
          std::size_t total = 0;
          for (const auto& ray : rays) {
            total += alternatives[k - 1]
                         ->candidates(ray.position, ray.direction)
                         .size();
          }
          candidates = static_cast<double>(total) / static_cast<double>(count);
        }
        std::cout << names[k] << " us/ray=" << timings[k][repeats / 2]
                  << " min=" << timings[k].front()
                  << " max=" << timings[k].back()
                  << " exact_tests/ray=" << candidates << '\n';
      }
    }
    // Check conservative fallback with unsupported and non-coaxial surfaces.
    owned.push_back(Surface::makeShared<PlaneSurface>(
        Transform3(Translation3(0, 0, 300)),
        std::make_shared<RectangleBounds>(100, 100)));
    owned.push_back(Surface::makeShared<CylinderSurface>(
        Transform3(Translation3(10, 20, 30)) *
            AngleAxis3(0.1, Vector3::UnitX()),
        std::make_shared<CylinderBounds>(50, 100)));
    surfaces = {owned[owned.size() - 2].get(), owned.back().get()};
    IntersectionMaterialAssigner fallbackBaseline(
        {surfaces, {}}, getDefaultLogger("fallback", Logging::ERROR));
    for (auto mode : {IndexedMaterialAssigner::Mode::Octree,
                      IndexedMaterialAssigner::Mode::RZ,
                      IndexedMaterialAssigner::Mode::LinearRZ,
                      IndexedMaterialAssigner::Mode::Angular}) {
      IndexedMaterialAssigner fallback(gctx, surfaces, mode);
      if (fallback.fallbackSize() != surfaces.size()) {
        throw std::runtime_error("Fallback eligibility failed");
      }
      compare(fallbackBaseline, fallback, gctx, edgeRays);
    }
    // Empty and mixed indexed/fallback lists exercise both traversal exits.
    IntersectionMaterialAssigner emptyBaseline(
        {{}, {}}, getDefaultLogger("empty", Logging::ERROR));
    surfaces.push_back(owned.front().get());
    IntersectionMaterialAssigner mixedBaseline(
        {surfaces, {}}, getDefaultLogger("mixed", Logging::ERROR));
    for (auto mode : {IndexedMaterialAssigner::Mode::Octree,
                      IndexedMaterialAssigner::Mode::RZ,
                      IndexedMaterialAssigner::Mode::LinearRZ,
                      IndexedMaterialAssigner::Mode::Angular}) {
      IndexedMaterialAssigner empty(gctx, {}, mode);
      compare(emptyBaseline, empty, gctx, edgeRays);
      IndexedMaterialAssigner mixed(gctx, surfaces, mode);
      compare(mixedBaseline, mixed, gctx, edgeRays);
    }
    cylinderExperiment(gctx, count, repeats);
    std::cout << "fallback_checks=passed checksum=" << checksum << '\n';
  } catch (const std::exception& e) {
    std::cerr << e.what() << '\n';
    return 1;
  }
}

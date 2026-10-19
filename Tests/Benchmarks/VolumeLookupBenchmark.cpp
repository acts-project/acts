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
#include <numeric>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

#include <Eigen/Geometry>

namespace {
using Clock = std::chrono::steady_clock;
using Volume = Acts::TrackingVolume;
using Vec = Acts::Vector3;
constexpr double tolerance = 1e-6;
constexpr std::size_t none = std::numeric_limits<std::size_t>::max();

// Benchmark-only broad phase. Exact Acts containment and sibling priority are
// retained. Build one median-split BVH for each sibling list.
struct Lookup {
  struct Node {
    Eigen::AlignedBox3d m_box;
    std::size_t m_begin;
    std::size_t m_end;
    std::size_t m_left = none;
    std::size_t m_right = none;
  };
  const Volume* m_volume;
  bool m_radial;
  bool m_octree;
  using OctBox = Acts::AxisAlignedBoundingBox<std::size_t, double, 3>;
  std::vector<std::unique_ptr<OctBox>> m_octStore;
  const OctBox* m_octRoot = nullptr;
  Eigen::AlignedBox3d m_box;
  std::vector<Lookup> m_children;
  std::vector<std::size_t> m_order;
  std::vector<Node> m_nodes;

  Lookup(const Lookup&) = delete;
  Lookup& operator=(const Lookup&) = delete;
  Lookup(Lookup&&) = default;
  Lookup& operator=(Lookup&&) = default;

  Lookup(const Acts::GeometryContext& gctx, const Volume& v, bool rz,
         bool oct = false)
      : m_volume(&v), m_radial(rz), m_octree(oct) {
    if (v.confinedVolumes() || !v.denseVolumes().empty()) {
      throw std::runtime_error("Prototype only supports Gen3 child volumes");
    }
    const auto bb = v.boundingBox(gctx, Vec::Constant(tolerance + 1e-9));
    m_box = Eigen::AlignedBox3d(bb.min(), bb.max());
    if (rz) {
      const auto* cylinder =
          dynamic_cast<const Acts::CylinderVolumeBounds*>(&v.volumeBounds());
      const auto& trf = v.localToGlobalTransform(gctx);
      if (cylinder == nullptr ||
          !trf.linear().isApprox(Acts::RotationMatrix3::Identity(), 1e-14) ||
          cylinder->get(Acts::CylinderVolumeBounds::eHalfPhiSector) !=
              std::numbers::pi) {
        throw std::runtime_error(
            "RZ prototype requires unrotated full cylinders");
      }
      // Offset cylinder axes are conservatively enclosed, including tiny
      // floating-point translations in the exported geometry.
      const double offset = trf.translation().head<2>().norm();
      m_box.min() = Vec(cylinder->get(Acts::CylinderVolumeBounds::eMinR) -
                            offset - tolerance - 1e-9,
                        bb.min().z(), -1);
      m_box.max() = Vec(cylinder->get(Acts::CylinderVolumeBounds::eMaxR) +
                            offset + tolerance + 1e-9,
                        bb.max().z(), 1);
    }
    for (const auto& child : v.volumes()) {
      m_children.emplace_back(gctx, child, rz, oct);
    }
    m_order.resize(m_children.size());
    std::iota(m_order.begin(), m_order.end(), 0);
    if (m_octree && !m_order.empty()) {
      std::vector<OctBox*> primitives;
      for (auto& i : m_order) {
        m_octStore.push_back(std::make_unique<OctBox>(
            &i, m_children[i].m_box.min(), m_children[i].m_box.max()));
        primitives.push_back(m_octStore.back().get());
      }
      m_octRoot =
          Acts::BoundingBoxHierarchy::makeOctree(m_octStore, primitives, 8);
    } else if (!m_order.empty()) {
      build(0, m_order.size());
    }
  }

  std::size_t build(std::size_t begin, std::size_t end) {
    const auto id = m_nodes.size();
    m_nodes.push_back({{}, begin, end});
    Eigen::AlignedBox3d centers;
    for (auto i = begin; i < end; ++i) {
      const auto& b = m_children[m_order[i]].m_box;
      m_nodes[id].m_box.extend(b);
      Vec c = (b.min() + b.max()) * 0.5;
      centers.extend(c);
    }
    if (end - begin > 4) {
      Eigen::Index axis = 0;
      (centers.max() - centers.min()).maxCoeff(&axis);
      const auto mid = begin + (end - begin) / 2;
      std::nth_element(m_order.begin() + begin, m_order.begin() + mid,
                       m_order.begin() + end, [&](auto a, auto b) {
                         return m_children[a].m_box.min()[axis] +
                                    m_children[a].m_box.max()[axis] <
                                m_children[b].m_box.min()[axis] +
                                    m_children[b].m_box.max()[axis];
                       });
      const auto left = build(begin, mid);
      const auto right = build(mid, end);
      m_nodes[id].m_left = left;
      m_nodes[id].m_right = right;
    }
    return id;
  }

  void search(const Acts::GeometryContext& gctx, std::size_t id, const Vec& p,
              const Vec& key, std::size_t& best) const {
    const auto& node = m_nodes[id];
    if (!node.m_box.contains(key)) {
      return;
    }
    if (node.m_left != none) {
      search(gctx, node.m_left, p, key, best);
      search(gctx, node.m_right, p, key, best);
      return;
    }
    for (auto i = node.m_begin; i < node.m_end; ++i) {
      const auto candidate = m_order[i];
      if (candidate < best && m_children[candidate].m_box.contains(key) &&
          m_children[candidate].m_volume->inside(gctx, p, tolerance)) {
        best = candidate;
      }
    }
  }

  const Volume* find(const Acts::GeometryContext& gctx, const Vec& p,
                     const Vec& key) const {
    if (!m_volume->inside(gctx, p, tolerance)) {
      return nullptr;
    }
    std::size_t best = none;
    if (m_octree) {
      Acts::BoundingBoxHierarchy::visitIntersecting(
          key, m_octRoot, [&](std::size_t candidate) {
            if (candidate < best &&
                m_children[candidate].m_volume->inside(gctx, p, tolerance)) {
              best = candidate;
            }
          });
    } else if (!m_nodes.empty()) {
      search(gctx, 0, p, key, best);
    }
    return best == none ? m_volume : m_children[best].find(gctx, p, key);
  }
  const Volume* operator()(const Acts::GeometryContext& gctx,
                           const Vec& p) const {
    return find(gctx, p, m_radial ? Vec(p.head<2>().norm(), p.z(), 0) : p);
  }
};

// Spatial partitions of the same conservative RZ rectangles as the BVH.
// Candidate lists retain sibling order; a rectangle may occupy multiple cells.
struct RzPartition {
  const Lookup& m_source;
  Eigen::AlignedBox3d m_domain;
  std::vector<RzPartition> m_children;
  std::array<std::vector<double>, 2> m_edges;
  std::vector<std::size_t> m_offsets;
  std::vector<std::size_t> m_candidates;

  explicit RzPartition(const Lookup& lookup) : m_source(lookup) {
    for (const auto& child : m_source.m_children) {
      m_children.emplace_back(child);
      m_domain.extend(child.m_box);
    }
    if (m_children.empty()) {
      return;
    }
    buildGrid();
  }

  std::size_t bin(int axis, double x) const {
    const auto& e = m_edges[axis];
    // Caller checks the domain. Include the final upper edge in the last bin.
    return std::min(
        static_cast<std::size_t>(std::upper_bound(e.begin(), e.end(), x) -
                                 e.begin() - 1),
        e.size() - 2);
  }

  void buildGrid() {
    for (int axis = 0; axis < 2; ++axis) {
      auto& e = m_edges[axis];
      for (const auto& child : m_source.m_children) {
        e.push_back(child.m_box.min()[axis]);
        e.push_back(child.m_box.max()[axis]);
      }
      std::sort(e.begin(), e.end());
      e.erase(std::unique(e.begin(), e.end()), e.end());
    }
    const auto nr = m_edges[0].size() - 1;
    const auto nz = m_edges[1].size() - 1;
    m_offsets.resize(nr * nz + 1);
    auto cells = [&](const Eigen::AlignedBox3d& b, auto&& visit) {
      for (auto r = bin(0, b.min()[0]); r <= bin(0, b.max()[0]); ++r) {
        for (auto z = bin(1, b.min()[1]); z <= bin(1, b.max()[1]); ++z) {
          visit(r * nz + z);
        }
      }
    };
    // Count, prefix sum, fill: contiguous candidate lists without per-bin
    // heaps.
    for (const auto& child : m_source.m_children) {
      cells(child.m_box, [&](auto cell) { ++m_offsets[cell + 1]; });
    }
    std::partial_sum(m_offsets.begin(), m_offsets.end(), m_offsets.begin());
    m_candidates.resize(m_offsets.back());
    auto cursor = m_offsets;
    for (std::size_t i = 0; i < m_children.size(); ++i) {
      cells(m_source.m_children[i].m_box,
            [&](auto cell) { m_candidates[cursor[cell]++] = i; });
    }
  }

  const Volume* find(const Acts::GeometryContext& gctx, const Vec& p,
                     const Vec& key) const {
    if (!m_source.m_volume->inside(gctx, p, tolerance)) {
      return nullptr;
    }
    if (!m_children.empty() && m_domain.contains(key)) {
      const auto cell =
          bin(0, key[0]) * (m_edges[1].size() - 1) + bin(1, key[1]);
      const auto begin = m_offsets[cell];
      const auto end = m_offsets[cell + 1];
      for (auto i = begin; i < end; ++i) {
        const auto candidate = m_candidates[i];
        const auto& child = m_source.m_children[candidate];
        if (child.m_box.contains(key) &&
            child.m_volume->inside(gctx, p, tolerance)) {
          return m_children[candidate].find(gctx, p, key);
        }
      }
    }
    return m_source.m_volume;
  }

  const Volume* operator()(const Acts::GeometryContext& gctx,
                           const Vec& p) const {
    return find(gctx, p, Vec(p.head<2>().norm(), p.z(), 0));
  }

  struct Stats {
    std::size_t m_cells = 0;
    std::size_t m_references = 0;
    std::size_t m_bytes = 0;
    std::size_t m_maxCandidates = 0;
  };
  void statistics(Stats& out) const {
    out.m_references += m_candidates.size();
    // Allocated array payload; excludes shared source boxes, object headers,
    // allocator overhead and temporary build storage.
    out.m_bytes +=
        sizeof(std::size_t) * (m_candidates.capacity() + m_offsets.capacity());
    for (const auto& e : m_edges) {
      out.m_bytes += sizeof(double) * e.capacity();
    }
    if (!m_offsets.empty()) {
      out.m_cells += m_offsets.size() - 1;
      for (std::size_t i = 1; i < m_offsets.size(); ++i) {
        out.m_maxCandidates =
            std::max(out.m_maxCandidates, m_offsets[i] - m_offsets[i - 1]);
      }
    }
    for (const auto& c : m_children) {
      c.statistics(out);
    }
  }
};

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
  std::unique_ptr<Lookup> xyz, rz, oct;
  std::cerr << "XYZ build ms: "
            << 1000 * elapsed([&] {
                 xyz = std::make_unique<Lookup>(gctx, world, false);
               })
            << '\n';
  std::cerr << "RZ build ms: "
            << 1000 * elapsed([&] {
                 rz = std::make_unique<Lookup>(gctx, world, true);
               })
            << '\n';
  std::cerr << "Octree build ms: "
            << 1000 * elapsed([&] {
                 oct = std::make_unique<Lookup>(gctx, world, false, true);
               })
            << '\n';
  std::unique_ptr<RzPartition> grid;
  std::cerr << "RZ grid build ms: "
            << 1000 *
                   elapsed([&] { grid = std::make_unique<RzPartition>(*rz); })
            << '\n';
  for (const auto* index : {grid.get()}) {
    RzPartition::Stats stats;
    index->statistics(stats);
    std::cerr << "RZ grid"
              << ": cells=" << stats.m_cells
              << ", references=" << stats.m_references
              << ", array_bytes=" << stats.m_bytes
              << ", max_candidates=" << stats.m_maxCandidates << '\n';
  }
  const auto gridLookup = [&](const Vec& p) { return (*grid)(gctx, p); };
  const auto octree = [&](const Vec& p) { return (*oct)(gctx, p); };
  const auto baseline = [&](const Vec& p) {
    return world.lowestTrackingVolume(gctx, p, tolerance);
  };
  const auto cartesian = [&](const Vec& p) { return (*xyz)(gctx, p); };
  const auto cylindrical = [&](const Vec& p) { return (*rz)(gctx, p); };
  auto productionGrid = converter.fromJson(gctx, payload);
  productionGrid->highestTrackingVolume()->setVolumeLookup(
      std::make_unique<Acts::RzGridVolumeLookup>(
          gctx, *productionGrid->highestTrackingVolume(),
          Acts::RzGridVolumeLookup::Config{}));
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
      if (gridLookup(p) != expected ||
          !sameVolume(productionGridLookup(p), expected)) {
        throw std::runtime_error("Lookup mismatch at partition edge");
      }
      ++edgeChecks;
    }
  };
  auto checkPartitions = [&](auto&& self, const RzPartition& index) -> void {
    for (double r : index.m_edges[0]) {
      for (double z : index.m_edges[1]) {
        checkEdge(r, z);
      }
    }
    for (const auto& child : index.m_children) {
      self(self, child);
    }
  };
  checkPartitions(checkPartitions, *grid);
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
      if (cartesian(p) != expected || cylindrical(p) != expected ||
          octree(p) != expected || gridLookup(p) != expected ||
          !sameVolume(productionGridLookup(p), expected)) {
        throw std::runtime_error("Lookup mismatch in " + workload);
      }
    }
    std::cerr << workload << ": validated " << points.size()
              << ", outside world " << misses << '\n';
    measure(workload, "linear", points, repeats, baseline);
    measure(workload, "xyz-octree", points, repeats, octree);
    measure(workload, "xyz-bvh", points, repeats, cartesian);
    measure(workload, "rz-bvh", points, repeats, cylindrical);
    measure(workload, "rz-grid", points, repeats, gridLookup);
    measure(workload, "production-rz-grid", points, repeats,
            productionGridLookup);
  }
  return 0;
} catch (const std::exception& e) {
  std::cerr << e.what() << '\n';
  return 1;
}

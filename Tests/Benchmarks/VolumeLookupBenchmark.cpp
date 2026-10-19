// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Geometry/CylinderVolumeBounds.hpp"
#include "Acts/Geometry/RzGridVolumeLookup.hpp"
#include "Acts/Geometry/RzKdTree4VolumeLookup.hpp"
#include "Acts/Geometry/RzKdTreeVolumeLookup.hpp"
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
    Eigen::AlignedBox3d box;
    std::size_t begin, end, left = none, right = none;
  };
  const Volume* volume;
  bool radial;
  bool octree;
  using OctBox = Acts::AxisAlignedBoundingBox<std::size_t, double, 3>;
  std::vector<std::unique_ptr<OctBox>> octStore;
  const OctBox* octRoot = nullptr;
  Eigen::AlignedBox3d box;
  std::vector<Lookup> children;
  std::vector<std::size_t> order;
  std::vector<Node> nodes;

  Lookup(const Lookup&) = delete;
  Lookup& operator=(const Lookup&) = delete;
  Lookup(Lookup&&) = default;
  Lookup& operator=(Lookup&&) = default;

  Lookup(const Acts::GeometryContext& gctx, const Volume& v, bool rz,
         bool oct = false)
      : volume(&v), radial(rz), octree(oct) {
    if (v.confinedVolumes() || !v.denseVolumes().empty()) {
      throw std::runtime_error("Prototype only supports Gen3 child volumes");
    }
    const auto bb = v.boundingBox(gctx, Vec::Constant(tolerance + 1e-9));
    box = Eigen::AlignedBox3d(bb.min(), bb.max());
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
      box.min() = Vec(cylinder->get(Acts::CylinderVolumeBounds::eMinR) -
                          offset - tolerance - 1e-9,
                      bb.min().z(), -1);
      box.max() = Vec(cylinder->get(Acts::CylinderVolumeBounds::eMaxR) +
                          offset + tolerance + 1e-9,
                      bb.max().z(), 1);
    }
    for (const auto& child : v.volumes()) {
      children.emplace_back(gctx, child, rz, oct);
    }
    order.resize(children.size());
    std::iota(order.begin(), order.end(), 0);
    if (octree && !order.empty()) {
      std::vector<OctBox*> primitives;
      for (auto& i : order) {
        octStore.push_back(std::make_unique<OctBox>(&i, children[i].box.min(),
                                                    children[i].box.max()));
        primitives.push_back(octStore.back().get());
      }
      octRoot = Acts::BoundingBoxHierarchy::makeOctree(octStore, primitives, 8);
    } else if (!order.empty()) {
      build(0, order.size());
    }
  }

  std::size_t build(std::size_t begin, std::size_t end) {
    const auto id = nodes.size();
    nodes.push_back({{}, begin, end});
    Eigen::AlignedBox3d centers;
    for (auto i = begin; i < end; ++i) {
      const auto& b = children[order[i]].box;
      nodes[id].box.extend(b);
      Vec c = (b.min() + b.max()) * 0.5;
      centers.extend(c);
    }
    if (end - begin > 4) {
      Eigen::Index axis = 0;
      (centers.max() - centers.min()).maxCoeff(&axis);
      const auto mid = begin + (end - begin) / 2;
      std::nth_element(
          order.begin() + begin, order.begin() + mid, order.begin() + end,
          [&](auto a, auto b) {
            return children[a].box.min()[axis] + children[a].box.max()[axis] <
                   children[b].box.min()[axis] + children[b].box.max()[axis];
          });
      const auto left = build(begin, mid);
      const auto right = build(mid, end);
      nodes[id].left = left;
      nodes[id].right = right;
    }
    return id;
  }

  void search(const Acts::GeometryContext& gctx, std::size_t id, const Vec& p,
              const Vec& key, std::size_t& best) const {
    const auto& node = nodes[id];
    if (!node.box.contains(key)) {
      return;
    }
    if (node.left != none) {
      search(gctx, node.left, p, key, best);
      search(gctx, node.right, p, key, best);
      return;
    }
    for (auto i = node.begin; i < node.end; ++i) {
      const auto candidate = order[i];
      if (candidate < best && children[candidate].box.contains(key) &&
          children[candidate].volume->inside(gctx, p, tolerance)) {
        best = candidate;
      }
    }
  }

  const Volume* find(const Acts::GeometryContext& gctx, const Vec& p,
                     const Vec& key) const {
    if (!volume->inside(gctx, p, tolerance)) {
      return nullptr;
    }
    std::size_t best = none;
    if (octree) {
      Acts::BoundingBoxHierarchy::visitIntersecting(
          key, octRoot, [&](std::size_t candidate) {
            if (candidate < best &&
                children[candidate].volume->inside(gctx, p, tolerance)) {
              best = candidate;
            }
          });
    } else if (!nodes.empty()) {
      search(gctx, 0, p, key, best);
    }
    return best == none ? volume : children[best].find(gctx, p, key);
  }
  const Volume* operator()(const Acts::GeometryContext& gctx,
                           const Vec& p) const {
    return find(gctx, p, radial ? Vec(p.head<2>().norm(), p.z(), 0) : p);
  }
};

// Spatial partitions of the same conservative RZ rectangles as the BVH.
// Candidate lists retain sibling order; a rectangle may occupy multiple cells.
struct RzPartition {
  struct Node {
    double split = 0;
    std::size_t left = none, right = none, begin = 0, end = 0;
    int axis = 0;
  };
  const Lookup& source;
  bool grid;
  Eigen::AlignedBox3d domain;
  std::vector<RzPartition> children;
  std::array<std::vector<double>, 2> edges;
  std::vector<std::size_t> offsets, candidates;
  std::vector<Node> nodes;

  RzPartition(const Lookup& lookup, bool useGrid)
      : source(lookup), grid(useGrid) {
    for (const auto& child : source.children) {
      children.emplace_back(child, grid);
      domain.extend(child.box);
    }
    if (children.empty()) {
      return;
    }
    if (grid) {
      buildGrid();
    } else {
      std::vector<std::size_t> ids(children.size());
      std::iota(ids.begin(), ids.end(), 0);
      buildKd(ids, 0);
    }
  }

  std::size_t bin(int axis, double x) const {
    const auto& e = edges[axis];
    // Caller checks the domain. Include the final upper edge in the last bin.
    return std::min(
        static_cast<std::size_t>(std::upper_bound(e.begin(), e.end(), x) -
                                 e.begin() - 1),
        e.size() - 2);
  }

  void buildGrid() {
    for (int axis = 0; axis < 2; ++axis) {
      auto& e = edges[axis];
      for (const auto& child : source.children) {
        e.push_back(child.box.min()[axis]);
        e.push_back(child.box.max()[axis]);
      }
      std::sort(e.begin(), e.end());
      e.erase(std::unique(e.begin(), e.end()), e.end());
    }
    const auto nr = edges[0].size() - 1;
    const auto nz = edges[1].size() - 1;
    offsets.resize(nr * nz + 1);
    auto cells = [&](const Eigen::AlignedBox3d& b, auto&& visit) {
      for (auto r = bin(0, b.min()[0]); r <= bin(0, b.max()[0]); ++r) {
        for (auto z = bin(1, b.min()[1]); z <= bin(1, b.max()[1]); ++z) {
          visit(r * nz + z);
        }
      }
    };
    // Count, prefix sum, fill: contiguous candidate lists without per-bin
    // heaps.
    for (const auto& child : source.children) {
      cells(child.box, [&](auto cell) { ++offsets[cell + 1]; });
    }
    std::partial_sum(offsets.begin(), offsets.end(), offsets.begin());
    candidates.resize(offsets.back());
    auto cursor = offsets;
    for (std::size_t i = 0; i < children.size(); ++i) {
      cells(source.children[i].box,
            [&](auto cell) { candidates[cursor[cell]++] = i; });
    }
  }

  std::size_t buildKd(const std::vector<std::size_t>& ids, std::size_t depth) {
    const auto index = nodes.size();
    nodes.emplace_back();
    int axis = -1;
    double split = 0;
    std::size_t cost = 2 * ids.size();
    if (ids.size() > 4 && depth < 24) {
      for (int a = 0; a < 2; ++a) {
        std::vector<double> centers;
        for (auto i : ids) {
          const auto& b = source.children[i].box;
          centers.push_back(std::midpoint(b.min()[a], b.max()[a]));
        }
        auto mid = centers.begin() + centers.size() / 2;
        std::nth_element(centers.begin(), mid, centers.end());
        const double cut = *mid;
        std::size_t nl = 0, nr = 0;
        for (auto i : ids) {
          const auto& b = source.children[i].box;
          if (b.min()[a] < cut) {
            ++nl;
          }
          if (b.max()[a] >= cut) {
            ++nr;
          }
        }
        // Both branches must reduce the candidate set; bound duplication and
        // stop degeneracy for coincident or strongly overlapping rectangles.
        if (nl > 0 && nr > 0 && nl < ids.size() && nr < ids.size() &&
            nl + nr < cost && nl + nr <= ids.size() + ids.size() / 2) {
          axis = a;
          split = cut;
          cost = nl + nr;
        }
      }
    }
    if (axis < 0) {
      nodes[index].begin = candidates.size();
      candidates.insert(candidates.end(), ids.begin(), ids.end());
      nodes[index].end = candidates.size();
    } else {
      std::vector<std::size_t> left, right;
      for (auto i : ids) {
        const auto& b = source.children[i].box;
        if (b.min()[axis] < split) {
          left.push_back(i);
        }
        if (b.max()[axis] >= split) {
          right.push_back(i);
        }
      }
      nodes[index].axis = axis;
      nodes[index].split = split;
      const auto l = buildKd(left, depth + 1);
      const auto r = buildKd(right, depth + 1);
      nodes[index].left = l;
      nodes[index].right = r;
    }
    return index;
  }

  const Volume* find(const Acts::GeometryContext& gctx, const Vec& p,
                     const Vec& key) const {
    if (!source.volume->inside(gctx, p, tolerance)) {
      return nullptr;
    }
    if (!children.empty() && domain.contains(key)) {
      std::size_t begin = 0;
      std::size_t end = 0;
      if (grid) {
        const auto cell =
            bin(0, key[0]) * (edges[1].size() - 1) + bin(1, key[1]);
        begin = offsets[cell];
        end = offsets[cell + 1];
      } else {
        std::size_t id = 0;
        while (nodes[id].left != none) {
          const auto& n = nodes[id];
          id = key[n.axis] < n.split ? n.left : n.right;
        }
        begin = nodes[id].begin;
        end = nodes[id].end;
      }
      for (auto i = begin; i < end; ++i) {
        const auto candidate = candidates[i];
        const auto& child = source.children[candidate];
        if (child.box.contains(key) &&
            child.volume->inside(gctx, p, tolerance)) {
          return children[candidate].find(gctx, p, key);
        }
      }
    }
    return source.volume;
  }

  const Volume* operator()(const Acts::GeometryContext& gctx,
                           const Vec& p) const {
    return find(gctx, p, Vec(p.head<2>().norm(), p.z(), 0));
  }

  struct Stats {
    std::size_t cells = 0, references = 0, bytes = 0, maxCandidates = 0;
  };
  void statistics(Stats& out) const {
    out.references += candidates.size();
    // Allocated array payload; excludes shared source boxes, object headers,
    // allocator overhead and temporary build storage.
    out.bytes +=
        sizeof(std::size_t) * (candidates.capacity() + offsets.capacity()) +
        sizeof(Node) * nodes.capacity();
    for (const auto& e : edges) {
      out.bytes += sizeof(double) * e.capacity();
    }
    if (grid && !offsets.empty()) {
      out.cells += offsets.size() - 1;
      for (std::size_t i = 1; i < offsets.size(); ++i) {
        out.maxCandidates =
            std::max(out.maxCandidates, offsets[i] - offsets[i - 1]);
      }
    } else {
      for (const auto& n : nodes) {
        if (n.left == none) {
          ++out.cells;
          out.maxCandidates = std::max(out.maxCandidates, n.end - n.begin);
        }
      }
    }
    for (const auto& c : children) {
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
  std::unique_ptr<RzPartition> kd, grid;
  std::cerr << "RZ kd build ms: "
            << 1000 * elapsed([&] {
                 kd = std::make_unique<RzPartition>(*rz, false);
               })
            << '\n';
  std::cerr << "RZ grid build ms: "
            << 1000 * elapsed([&] {
                 grid = std::make_unique<RzPartition>(*rz, true);
               })
            << '\n';
  for (const auto* index : {kd.get(), grid.get()}) {
    RzPartition::Stats stats;
    index->statistics(stats);
    std::cerr << (index->grid ? "RZ grid" : "RZ kd")
              << ": cells=" << stats.cells
              << ", references=" << stats.references
              << ", array_bytes=" << stats.bytes
              << ", max_candidates=" << stats.maxCandidates << '\n';
  }
  const auto kdLookup = [&](const Vec& p) { return (*kd)(gctx, p); };
  const auto gridLookup = [&](const Vec& p) { return (*grid)(gctx, p); };
  const auto octree = [&](const Vec& p) { return (*oct)(gctx, p); };
  const auto baseline = [&](const Vec& p) {
    return world.lowestTrackingVolume(gctx, p, tolerance);
  };
  const auto cartesian = [&](const Vec& p) { return (*xyz)(gctx, p); };
  const auto cylindrical = [&](const Vec& p) { return (*rz)(gctx, p); };
  auto productionKd = converter.fromJson(gctx, payload);
  auto productionGrid = converter.fromJson(gctx, payload);
  auto existingKdGeometry = converter.fromJson(gctx, payload);
  std::cerr << "Existing 4D KD build ms: "
            << 1000 * elapsed([&] {
                 auto& root = *existingKdGeometry->highestTrackingVolume();
                 root.setVolumeLookup(
                     std::make_unique<Acts::RzKdTree4VolumeLookup>(
                         gctx, root, Acts::RzKdTree4VolumeLookup::Config{}));
               })
            << '\n';
  const auto existingKdLookup = [&](const Vec& p) {
    return *existingKdGeometry->resolveLowestTrackingVolume(
        gctx, p, std::nullopt, nullptr, tolerance);
  };
  productionKd->highestTrackingVolume()->setVolumeLookup(
      std::make_unique<Acts::RzKdTreeVolumeLookup>(
          gctx, *productionKd->highestTrackingVolume(),
          Acts::RzKdTreeVolumeLookup::Config{}));
  productionGrid->highestTrackingVolume()->setVolumeLookup(
      std::make_unique<Acts::RzGridVolumeLookup>(
          gctx, *productionGrid->highestTrackingVolume(),
          Acts::RzGridVolumeLookup::Config{}));
  const auto sameVolume = [](const Volume* lhs, const Volume* rhs) {
    return lhs == nullptr || rhs == nullptr
               ? lhs == rhs
               : lhs->geometryId() == rhs->geometryId();
  };
  const auto productionKdLookup = [&](const Vec& p) {
    return *productionKd->resolveLowestTrackingVolume(gctx, p, std::nullopt,
                                                      nullptr, tolerance);
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
      if (kdLookup(p) != expected || gridLookup(p) != expected ||
          !sameVolume(productionKdLookup(p), expected) ||
          !sameVolume(productionGridLookup(p), expected) ||
          !sameVolume(existingKdLookup(p), expected)) {
        throw std::runtime_error("Lookup mismatch at partition edge");
      }
      ++edgeChecks;
    }
  };
  auto checkPartitions = [&](auto&& self, const RzPartition& index) -> void {
    if (index.grid) {
      for (double r : index.edges[0]) {
        for (double z : index.edges[1]) {
          checkEdge(r, z);
        }
      }
    } else {
      for (const auto& node : index.nodes) {
        if (node.left == none) {
          continue;
        }
        for (const auto& child : index.source.children) {
          Vec key = (child.box.min() + child.box.max()) * 0.5;
          key[node.axis] = node.split;
          checkEdge(key[0], key[1]);
        }
      }
    }
    for (const auto& child : index.children) {
      self(self, child);
    }
  };
  checkPartitions(checkPartitions, *grid);
  checkPartitions(checkPartitions, *kd);
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
          octree(p) != expected || kdLookup(p) != expected ||
          gridLookup(p) != expected ||
          !sameVolume(productionKdLookup(p), expected) ||
          !sameVolume(productionGridLookup(p), expected) ||
          !sameVolume(existingKdLookup(p), expected)) {
        throw std::runtime_error("Lookup mismatch in " + workload);
      }
    }
    std::cerr << workload << ": validated " << points.size()
              << ", outside world " << misses << '\n';
    measure(workload, "linear", points, repeats, baseline);
    measure(workload, "xyz-octree", points, repeats, octree);
    measure(workload, "xyz-bvh", points, repeats, cartesian);
    measure(workload, "rz-bvh", points, repeats, cylindrical);
    measure(workload, "rz-kd", points, repeats, kdLookup);
    measure(workload, "rz-grid", points, repeats, gridLookup);
    measure(workload, "existing-4d-kd", points, repeats, existingKdLookup);
    measure(workload, "production-rz-kd", points, repeats, productionKdLookup);
    measure(workload, "production-rz-grid", points, repeats,
            productionGridLookup);
  }
  return 0;
} catch (const std::exception& e) {
  std::cerr << e.what() << '\n';
  return 1;
}

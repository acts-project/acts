// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Material/interface/IAssignmentFinder.hpp"
#include "Acts/Surfaces/CylinderBounds.hpp"
#include "Acts/Surfaces/CylinderSurface.hpp"
#include "Acts/Surfaces/DiscSurface.hpp"
#include "Acts/Surfaces/RadialBounds.hpp"
#include "Acts/Utilities/BoundingBox.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <memory>
#include <numeric>
#include <vector>

#include "OnAxisCylinderIntersection.hpp"

namespace ActsTests {

/// Benchmark-only surface assigner. Bounds are valid for the construction
/// geometry context only: rebuild after alignment changes. Unsupported surfaces
/// are always intersected. Does not provide volume assignments.
class IndexedMaterialAssigner final : public Acts::IAssignmentFinder {
 public:
  enum class Mode { Octree, RZ, LinearRZ, Angular, Brute };

  IndexedMaterialAssigner(const Acts::GeometryContext& gctx,
                          std::vector<const Acts::Surface*> surfaces, Mode mode,
                          bool fastCylinder = false)
      : m_surfaces(std::move(surfaces)),
        m_mode(mode),
        m_fastCylinder(fastCylinder) {
    m_bounds.resize(m_surfaces.size());
    for (const auto* surface : m_surfaces) {
      m_cylinders.push_back(
          dynamic_cast<const Acts::CylinderSurface*>(surface));
    }
    std::vector<Box*> primitives;
    std::vector<std::size_t> indices;
    for (std::size_t i = 0; i < m_surfaces.size(); ++i) {
      const auto& surface = *m_surfaces[i];
      const auto& trf = surface.localToGlobalTransform(gctx);
      // Exact eligibility check: a near-coaxial surface is not necessarily
      // enclosed by these bounds. Arbitrary in-plane rotation is fine.
      if (trf.translation().x() != 0 || trf.translation().y() != 0 ||
          trf.linear()(0, 2) != 0 || trf.linear()(1, 2) != 0 ||
          std::abs(trf.linear()(2, 2)) != 1) {
        m_fallback.push_back(i);
        continue;
      }
      Bounds b{};
      const double z = trf.translation().z();
      if (auto* cylinder =
              dynamic_cast<const Acts::CylinderSurface*>(&surface)) {
        const auto& cb = cylinder->bounds();
        b = {cb.get(Acts::CylinderBounds::eR), cb.get(Acts::CylinderBounds::eR),
             z - cb.get(Acts::CylinderBounds::eHalfLengthZ),
             z + cb.get(Acts::CylinderBounds::eHalfLengthZ)};
      } else if (dynamic_cast<const Acts::DiscSurface*>(&surface) != nullptr) {
        auto* rb = dynamic_cast<const Acts::RadialBounds*>(&surface.bounds());
        if (rb == nullptr) {
          m_fallback.push_back(i);
          continue;
        }
        b = {rb->get(Acts::RadialBounds::eMinR),
             rb->get(Acts::RadialBounds::eMaxR), z, z};
      } else {
        m_fallback.push_back(i);
        continue;
      }
      // Match/exceed the exact surface intersection's boundary tolerance and
      // make planar/cylindrical shells non-degenerate for conservative pruning.
      constexpr double padding = 2 * Acts::s_onSurfaceTolerance;
      b.rmin = std::max(0., b.rmin - padding);
      b.rmax += padding;
      b.zmin -= padding;
      b.zmax += padding;
      m_bounds[i] = b;
      indices.push_back(i);
      if (mode == Mode::Octree) {
        m_boxes.push_back(std::make_unique<Box>(
            &m_surfaces[i], Acts::Vector3(-b.rmax, -b.rmax, b.zmin),
            Acts::Vector3(b.rmax, b.rmax, b.zmax)));
        primitives.push_back(m_boxes.back().get());
      }
    }
    m_indexed = indices;
    if (mode == Mode::Angular) {
      // Fixed origin (0,0,0): each conservative r-z rectangle projects to
      // a cot(theta) interval. Keep singular radial bounds unconditionally.
      std::vector<std::pair<double, double>> intervals;
      for (auto i : indices) {
        const auto& b = m_bounds[i];
        if (b.rmin == 0) {
          m_angularAlways.push_back(i);
          intervals.push_back({0, 0});
          continue;
        }
        const std::array<double, 4> slopes{b.zmin / b.rmin, b.zmin / b.rmax,
                                           b.zmax / b.rmin, b.zmax / b.rmax};
        auto [lo, hi] = std::minmax_element(slopes.begin(), slopes.end());
        intervals.push_back({*lo, *hi});
        m_endpoints.push_back(*lo);
        m_endpoints.push_back(*hi);
      }
      std::ranges::sort(m_endpoints);
      m_endpoints.erase(std::unique(m_endpoints.begin(), m_endpoints.end()),
                        m_endpoints.end());
      m_angularCells.resize(m_endpoints.size() + 1);
      for (std::size_t k = 0; k < indices.size(); ++k) {
        if (m_bounds[indices[k]].rmin == 0) {
          continue;
        }
        auto [lo, hi] = intervals[k];
        auto begin =
            std::lower_bound(m_endpoints.begin(), m_endpoints.end(), lo) -
            m_endpoints.begin() + 1;
        auto end =
            std::lower_bound(m_endpoints.begin(), m_endpoints.end(), hi) -
            m_endpoints.begin() + 1;
        for (auto cell = begin; cell < end; ++cell) {
          m_angularCells[cell].push_back(indices[k]);
        }
      }
    }
    if (mode == Mode::Octree && !primitives.empty()) {
      m_root = Acts::BoundingBoxHierarchy::makeOctree(m_boxes, primitives, 6);
    } else if (mode == Mode::RZ && !indices.empty()) {
      build(indices, 0, indices.size());
    }
  }

  std::vector<std::size_t> candidates(const Acts::Vector3& p,
                                      const Acts::Vector3& d) const {
    if (m_mode == Mode::Brute) {
      std::vector<std::size_t> all(m_surfaces.size());
      std::iota(all.begin(), all.end(), 0);
      return all;
    }
    std::vector<std::size_t> result = m_fallback;
    if (m_mode == Mode::Angular && p.isZero(0) &&
        d.head<2>().squaredNorm() > 0) {
      const double slope = d.z() / d.head<2>().norm();
      auto it = std::lower_bound(m_endpoints.begin(), m_endpoints.end(), slope);
      const auto cell = static_cast<std::size_t>(it - m_endpoints.begin());
      result.insert(result.end(), m_angularAlways.begin(),
                    m_angularAlways.end());
      result.insert(result.end(), m_angularCells[cell].begin(),
                    m_angularCells[cell].end());
      // At an endpoint include both neighbours, retaining closed intervals.
      if (it != m_endpoints.end() && *it == slope) {
        result.insert(result.end(), m_angularCells[cell + 1].begin(),
                      m_angularCells[cell + 1].end());
      }
      std::ranges::sort(result);
      result.erase(std::unique(result.begin(), result.end()), result.end());
      return result;
    }
    if (m_mode == Mode::Octree) {
      // A robust inclusive slab predicate avoids 0*infinity on box boundaries.
      Acts::Vector3 inverse;
      for (int axis = 0; axis < 3; ++axis) {
        inverse[axis] = d[axis] == 0 ? 0 : 1 / d[axis];
      }
      const Box* node = m_root;
      while (node != nullptr) {
        double lo = 0, hi = std::numeric_limits<double>::infinity();
        bool hit = true;
        for (int axis = 0; axis < 3; ++axis) {
          if (d[axis] == 0) {
            hit &= p[axis] >= node->min()[axis] && p[axis] <= node->max()[axis];
          } else {
            double a = (node->min()[axis] - p[axis]) * inverse[axis];
            double b = (node->max()[axis] - p[axis]) * inverse[axis];
            lo = std::max(lo, std::min(a, b));
            hi = std::min(hi, std::max(a, b));
          }
        }
        if (!hit || lo > hi) {
          node = node->getSkip();
        } else if (!node->hasEntity()) {
          node = node->getLeftChild();
        } else {
          result.push_back(node->entity() - m_surfaces.data());
          node = node->getSkip();
        }
      }
    } else {
      const RadialRay ray(p, d);
      if (m_mode == Mode::LinearRZ || m_mode == Mode::Angular) {
        for (auto i : m_indexed) {
          if (ray.intersects(m_bounds[i])) {
            result.push_back(i);
          }
        }
      } else {
        std::size_t i = 0;
        while (i < m_nodes.size()) {
          const auto& node = m_nodes[i];
          if (!ray.intersects(node.bounds)) {
            i = node.skip;
          } else {
            if (node.surface != invalid) {
              result.push_back(node.surface);
            }
            ++i;
          }
        }
      }
    }
    // Retain original input order, including tie behaviour of the exact
    // assigner's subsequent path-length sort.
    std::ranges::sort(result);
    return result;
  }

  std::pair<std::vector<SurfaceAssignment>, std::vector<VolumeAssignment>>
  assignmentCandidates(const Acts::GeometryContext& gctx,
                       const Acts::MagneticFieldContext& /*mctx*/,
                       const Acts::Vector3& p,
                       const Acts::Vector3& d) const override {
    struct Hit {
      Acts::Intersection3D intersection;
      const Acts::Surface* surface;
    };
    std::vector<Hit> hits;
    for (auto i : candidates(p, d)) {
      const auto* cylinder = m_fastCylinder ? m_cylinders[i] : nullptr;
      auto multi =
          cylinder != nullptr
              ? intersectOnAxisCylinder(*cylinder, gctx, p, d,
                                        Acts::BoundaryTolerance::None())
              : m_surfaces[i]->intersect(gctx, p, d,
                                         Acts::BoundaryTolerance::None());
      const auto& hit = multi.closestForward();
      if (hit.status() >= Acts::IntersectionStatus::reachable &&
          hit.pathLength() > 0) {
        hits.push_back({hit, m_surfaces[i]});
      }
    }
    std::ranges::sort(hits, [](const Hit& a, const Hit& b) {
      return Acts::Intersection3D::pathLengthOrder(a.intersection,
                                                   b.intersection);
    });
    std::vector<SurfaceAssignment> result;
    result.reserve(hits.size());
    for (const auto& hit : hits) {
      result.push_back({hit.surface, hit.intersection.position(), d});
    }
    return {std::move(result), {}};
  }

  std::size_t fallbackSize() const { return m_fallback.size(); }

 private:
  struct Bounds {
    double rmin, rmax, zmin, zmax;
  };
  // Over a z-clipped forward ray interval, r(t)^2 is a convex quadratic.
  // Its minimum is at the closest transverse approach or an interval endpoint;
  // its maximum is at an endpoint (or infinity). Overlap with [rmin,rmax]
  // proves that the ray can reach the annular cylindrical bounding volume.
  struct RadialRay {
    Acts::Vector3 p, d;
    double a, closest, inverseZ, inverseR;
    bool onAxis;
    RadialRay(const Acts::Vector3& origin, const Acts::Vector3& direction)
        : p(origin),
          d(direction),
          a(d.head<2>().squaredNorm()),
          closest(a == 0 ? 0 : -p.head<2>().dot(d.head<2>()) / a),
          inverseZ(d.z() == 0 ? 0 : 1 / d.z()),
          inverseR(a == 0 ? 0 : 1 / std::sqrt(a)),
          onAxis(p.x() == 0 && p.y() == 0) {}
    bool intersects(const Bounds& b) const {
      double lo = 0, hi = std::numeric_limits<double>::infinity();
      if (d.z() == 0) {
        if (p.z() < b.zmin || p.z() > b.zmax) {
          return false;
        }
      } else {
        double t0 = (b.zmin - p.z()) * inverseZ;
        double t1 = (b.zmax - p.z()) * inverseZ;
        lo = std::max(0., std::min(t0, t1));
        hi = std::max(t0, t1);
        if (lo > hi) {
          return false;
        }
      }
      // On the symmetry axis, r(t) = t * |d_xy| for every forward ray:
      // the curved 3D query becomes an ordinary 2D r-z slab test.
      if (onAxis && a != 0) {
        return std::max(lo, b.rmin * inverseR) <=
               std::min(hi, b.rmax * inverseR);
      }
      auto radius2 = [&](double t) {
        return (p.head<2>() + t * d.head<2>()).squaredNorm();
      };
      if (a == 0) {
        const double r2 = p.head<2>().squaredNorm();
        return r2 >= b.rmin * b.rmin && r2 <= b.rmax * b.rmax;
      }
      const double rlo = radius2(std::clamp(closest, lo, hi));
      const double rhi = std::isfinite(hi)
                             ? std::max(radius2(lo), radius2(hi))
                             : std::numeric_limits<double>::infinity();
      return rlo <= b.rmax * b.rmax && rhi >= b.rmin * b.rmin;
    }
  };
  static constexpr std::size_t invalid =
      std::numeric_limits<std::size_t>::max();
  struct Node {
    Bounds bounds;
    std::size_t skip, surface = invalid;
  };
  void build(std::vector<std::size_t>& indices, std::size_t begin,
             std::size_t end) {
    Bounds b = m_bounds[indices[begin]];
    for (std::size_t i = begin + 1; i < end; ++i) {
      const auto& other = m_bounds[indices[i]];
      b = {std::min(b.rmin, other.rmin), std::max(b.rmax, other.rmax),
           std::min(b.zmin, other.zmin), std::max(b.zmax, other.zmax)};
    }
    const auto node = m_nodes.size();
    m_nodes.push_back({b, 0});
    if (end - begin == 1) {
      m_nodes[node].surface = indices[begin];
    } else {
      const bool splitZ = b.zmax - b.zmin > b.rmax - b.rmin;
      const auto mid = begin + (end - begin) / 2;
      std::nth_element(indices.begin() + begin, indices.begin() + mid,
                       indices.begin() + end, [&](auto i, auto j) {
                         const auto& x = m_bounds[i];
                         const auto& y = m_bounds[j];
                         return splitZ ? x.zmin + x.zmax < y.zmin + y.zmax
                                       : x.rmin + x.rmax < y.rmin + y.rmax;
                       });
      build(indices, begin, mid);
      build(indices, mid, end);
    }
    m_nodes[node].skip = m_nodes.size();
  }
  using Box = Acts::AxisAlignedBoundingBox<const Acts::Surface*, double, 3>;
  std::vector<const Acts::Surface*> m_surfaces;
  std::vector<const Acts::CylinderSurface*> m_cylinders;
  Mode m_mode;
  bool m_fastCylinder;
  std::vector<double> m_endpoints;
  std::vector<std::vector<std::size_t>> m_angularCells;
  std::vector<std::size_t> m_angularAlways;
  std::vector<Bounds> m_bounds;
  std::vector<std::size_t> m_fallback, m_indexed;
  std::vector<Node> m_nodes;
  std::vector<std::unique_ptr<Box>> m_boxes;
  const Box* m_root = nullptr;
};
}  // namespace ActsTests

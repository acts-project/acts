// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Geometry/RzKdTreeVolumeLookup.hpp"

#include <algorithm>
#include <array>
#include <numeric>
#include <span>
#include <stdexcept>

#include "VolumeLookupImpl.hpp"

namespace Acts {
namespace {
// Acts::KDTree stores points and supports point range searches. This partition
// instead stores overlapping rectangles in every intersected branch, so a point
// query follows one branch and returns candidates in hierarchy priority order.
struct Partition {
  using Config = RzKdTreeVolumeLookup::Config;
  struct Node {
    double split = 0;
    std::size_t left = detail::noRzNode;
    std::size_t right = detail::noRzNode;
    std::size_t begin = 0;
    std::size_t end = 0;
    int axis = 0;
  };
  std::vector<Node> nodes;
  std::vector<std::size_t> m_candidates;
  template <typename Children>
  bool build(const Children& children, const Config& config) {
    std::vector<std::size_t> ids;
    for (std::size_t i = 0; i < children.size(); ++i) {
      if (children[i].indexable) {
        ids.push_back(i);
      }
    }
    buildKd(children, config, ids, 0);
    return true;
  }
  template <typename Children>
  std::size_t buildKd(const Children& children, const Config& config,
                      const std::vector<std::size_t>& ids, std::size_t depth) {
    const auto index = nodes.size();
    nodes.emplace_back();
    int axis = -1;
    double split = 0;
    std::size_t cost = 2 * ids.size();
    if (ids.size() > config.leafSize && depth < config.maxDepth) {
      for (int a = 0; a < 2; ++a) {
        std::vector<double> centers;
        for (auto i : ids) {
          const auto& b = children[i].box;
          centers.push_back(std::midpoint(b.min()[a], b.max()[a]));
        }
        auto mid = centers.begin() + centers.size() / 2;
        std::nth_element(centers.begin(), mid, centers.end());
        const double cut = *mid;
        std::size_t nl = 0;
        std::size_t nr = 0;
        for (auto i : ids) {
          const auto& b = children[i].box;
          nl += b.min()[a] < cut;
          nr += b.max()[a] >= cut;
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
      nodes[index].begin = m_candidates.size();
      m_candidates.insert(m_candidates.end(), ids.begin(), ids.end());
      nodes[index].end = m_candidates.size();
    } else {
      std::vector<std::size_t> left;
      std::vector<std::size_t> right;
      for (auto i : ids) {
        const auto& b = children[i].box;
        if (b.min()[axis] < split) {
          left.push_back(i);
        }
        if (b.max()[axis] >= split) {
          right.push_back(i);
        }
      }
      nodes[index].axis = axis;
      nodes[index].split = split;
      const auto l = buildKd(children, config, left, depth + 1);
      const auto r = buildKd(children, config, right, depth + 1);
      nodes[index].left = l;
      nodes[index].right = r;
    }
    return index;
  }

  template <typename Visitor>
  void visitCandidates(const Vector2& key, Visitor&& visitor) const {
    for (auto id : candidates(key)) {
      if (!visitor(id)) {
        break;
      }
    }
  }

  std::span<const std::size_t> candidates(const Vector2& key) const {
    std::size_t id = 0;
    while (nodes[id].left != detail::noRzNode) {
      const auto& node = nodes[id];
      id = key[node.axis] < node.split ? node.left : node.right;
    }
    return std::span{m_candidates}.subspan(nodes[id].begin,
                                           nodes[id].end - nodes[id].begin);
  }
};
}  // namespace
struct RzKdTreeVolumeLookup::Impl : detail::RzLookupIndex<Partition> {
  Config config;
  Impl(const GeometryContext& gctx, const TrackingVolume& volume,
       const Config& cfg)
      : detail::RzLookupIndex<Partition>(gctx, volume, cfg), config(cfg) {}
};

RzKdTreeVolumeLookup::RzKdTreeVolumeLookup(const GeometryContext& gctx,
                                           const TrackingVolume& world,
                                           const Config& config,
                                           const Logger& logger) {
  if (!std::isfinite(config.fallbackWarningFraction) ||
      config.fallbackWarningFraction < 0 ||
      config.fallbackWarningFraction > 1 ||
      !std::isfinite(config.maxTolerance) || config.maxTolerance < 0 ||
      config.leafSize == 0 || config.maxDepth == 0 || config.maxDepth > 64) {
    throw std::invalid_argument("Invalid RZ volume lookup configuration");
  }
  m_impl = std::make_unique<Impl>(gctx, world, config);
  m_impl->warnAboutFallback(world, config.fallbackWarningFraction, logger);
}

RzKdTreeVolumeLookup::~RzKdTreeVolumeLookup() = default;

const RzKdTreeVolumeLookup::Config& RzKdTreeVolumeLookup::config() const {
  return m_impl->config;
}

const TrackingVolume* RzKdTreeVolumeLookup::lowestTrackingVolume(
    const TrackingVolume& root, const GeometryContext& gctx,
    const Vector3& position, double tolerance) const {
  if (!std::isfinite(tolerance) || tolerance < 0 ||
      tolerance > m_impl->config.maxTolerance) {
    return detail::linearVolumeScan(root, gctx, position, tolerance);
  }
  return m_impl->find(root, gctx, position,
                      Vector2(position.head<2>().norm(), position.z()),
                      tolerance);
}
}  // namespace Acts

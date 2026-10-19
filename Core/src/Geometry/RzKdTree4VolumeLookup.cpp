// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Geometry/RzKdTree4VolumeLookup.hpp"

#include "Acts/Utilities/KDTree.hpp"
#include "Acts/Utilities/VectorHelpers.hpp"

#include <limits>
#include <stdexcept>

#include "VolumeLookupImpl.hpp"

namespace Acts {
namespace {
// Rectangle containment is an orthogonal range query on four endpoint keys.
struct Partition {
  using Config = RzKdTree4VolumeLookup::Config;
  using Tree = KDTree<4, std::size_t>;
  std::unique_ptr<Tree> tree;

  template <typename Children>
  bool build(const Children& children, const Config& /*config*/) {
    Tree::vector_t points;
    for (std::size_t i = 0; i < children.size(); ++i) {
      if (children[i].indexable) {
        const auto& box = children[i].box;
        points.push_back(
            {{box.min()[0], box.min()[1], box.max()[0], box.max()[1]}, i});
      }
    }
    tree = std::make_unique<Tree>(std::move(points));
    return true;
  }

  template <typename Visitor>
  void visitCandidates(const Vector2& key, Visitor&& visitor) const {
    const double infinity = std::numeric_limits<double>::infinity();
    // RangeXD is half-open; include equality for the minimum endpoints.
    Tree::range_t range({-infinity, -infinity, key[0], key[1]},
                        {std::nextafter(key[0], infinity),
                         std::nextafter(key[1], infinity), infinity, infinity});
    tree->rangeSearchMapDiscard(range, [&](const auto&, std::size_t id) {
      // Tree traversal is unordered, so a match cannot terminate the search.
      visitor(id);
    });
  }
};
}  // namespace
struct RzKdTree4VolumeLookup::Impl : detail::RzLookupIndex<Partition> {
  Config config;
  Impl(const GeometryContext& gctx, const TrackingVolume& volume,
       const Config& cfg)
      : detail::RzLookupIndex<Partition>(gctx, volume, cfg), config(cfg) {}
};

RzKdTree4VolumeLookup::RzKdTree4VolumeLookup(const GeometryContext& gctx,
                                             const TrackingVolume& world,
                                             const Config& config,
                                             const Logger& logger) {
  if (!std::isfinite(config.fallbackWarningFraction) ||
      config.fallbackWarningFraction < 0 ||
      config.fallbackWarningFraction > 1 ||
      !std::isfinite(config.maxTolerance) || config.maxTolerance < 0) {
    throw std::invalid_argument("Invalid RZ volume lookup configuration");
  }
  m_impl = std::make_unique<Impl>(gctx, world, config);
  m_impl->warnAboutFallback(world, config.fallbackWarningFraction, logger);
}

RzKdTree4VolumeLookup::~RzKdTree4VolumeLookup() = default;

const RzKdTree4VolumeLookup::Config& RzKdTree4VolumeLookup::config() const {
  return m_impl->config;
}

const TrackingVolume* RzKdTree4VolumeLookup::lowestTrackingVolume(
    const TrackingVolume& root, const GeometryContext& gctx,
    const Vector3& position, double tolerance) const {
  if (!std::isfinite(tolerance) || tolerance < 0 ||
      tolerance > m_impl->config.maxTolerance) {
    return detail::linearVolumeScan(root, gctx, position, tolerance);
  }
  return m_impl->find(root, gctx, position,
                      Vector2(VectorHelpers::perp(position), position.z()),
                      tolerance);
}
}  // namespace Acts

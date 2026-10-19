// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Geometry/RzGridVolumeLookup.hpp"

#include "Acts/Utilities/Axis.hpp"
#include "Acts/Utilities/MultiAxis.hpp"
#include "Acts/Utilities/VectorHelpers.hpp"

#include <algorithm>
#include <array>
#include <numeric>
#include <optional>
#include <span>
#include <stdexcept>

#include "VolumeLookupImpl.hpp"

namespace Acts {
namespace {
// Nonuniform grid: collect all rectangle edges separately in r and z, sort
// and deduplicate them, then form their Cartesian product. A rectangle is
// assigned to every bin it overlaps, including the bin at its upper edge so
// a query exactly on that edge still finds it.
//
// Candidate lists use compressed sparse row storage: offsets[cell] and
// offsets[cell + 1] delimit a cell's slice of m_candidates. Count assignments,
// prefix-sum the counts, then fill slices in hierarchy order. Queries use two
// binary searches and scan one slice; shared code performs exact containment.
struct Partition {
  using Config = RzGridVolumeLookup::Config;
  using RzAxis = Axis<AxisType::Variable, AxisBoundaryType::Bound>;
  using RzAxes = MultiAxis<RzAxis, RzAxis>;
  std::optional<RzAxes> axes;
  std::vector<std::size_t> offsets;
  std::vector<std::size_t> m_candidates;

  template <typename Children>
  bool build(const Children& children, const Config& config) {
    std::array<std::vector<double>, 2> edges;
    for (int axis = 0; axis < 2; ++axis) {
      auto& e = edges[axis];
      for (const auto& child : children) {
        if (!child.indexable) {
          continue;
        }
        e.push_back(child.box.min()[axis]);
        e.push_back(child.box.max()[axis]);
      }
      std::sort(e.begin(), e.end());
      e.erase(std::unique(e.begin(), e.end()), e.end());
    }
    // MultiAxis global indices reserve under-/overflow slots on each axis,
    // even with bounded axes. Include these slots in the allocation limit.
    const auto nr = edges[0].size() + 1;
    const auto nz = edges[1].size() + 1;
    if (nr > config.maxGridCells / nz) {
      return false;
    }
    axes.emplace(RzAxis(std::move(edges[0])), RzAxis(std::move(edges[1])));
    offsets.resize(axes->getNTotalBins(true) + 1);
    auto cells = [&](const Eigen::AlignedBox2d& b, auto&& visit) {
      const auto lower = axes->getLocalBinsFromPoint(b.min());
      const auto upper = axes->getLocalBinsFromPoint(b.max());
      // Include the bin at the upper edge. Bound axes clamp the final edge to
      // the last regular bin, matching the inclusive rectangle containment.
      for (auto r = lower[0]; r <= upper[0]; ++r) {
        for (auto z = lower[1]; z <= upper[1]; ++z) {
          visit(axes->getGlobalBinFromLocalBins({r, z}));
        }
      }
    };
    // Count, prefix sum, fill: contiguous candidate lists without per-bin
    // heaps.
    for (const auto& child : children) {
      if (!child.indexable) {
        continue;
      }
      cells(child.box, [&](auto cell) { ++offsets[cell + 1]; });
    }
    std::partial_sum(offsets.begin(), offsets.end(), offsets.begin());
    m_candidates.resize(offsets.back());
    auto cursor = offsets;
    for (std::size_t i = 0; i < children.size(); ++i) {
      if (!children[i].indexable) {
        continue;
      }
      cells(children[i].box,
            [&](auto cell) { m_candidates[cursor[cell]++] = i; });
    }
    return true;
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
    // The caller checks the inclusive domain; MultiAxis::isInside() would
    // exclude the upper edge, while Bound axes alone also clamp outside points.
    const auto cell = axes->getGlobalBinFromPoint(key);
    return std::span{m_candidates}.subspan(offsets[cell],
                                           offsets[cell + 1] - offsets[cell]);
  }
};
}  // namespace
struct RzGridVolumeLookup::Impl : detail::RzLookupIndex<Partition> {
  Config config;
  Impl(const GeometryContext& gctx, const TrackingVolume& volume,
       const Config& cfg)
      : detail::RzLookupIndex<Partition>(gctx, volume, cfg), config(cfg) {}
};

RzGridVolumeLookup::RzGridVolumeLookup(const GeometryContext& gctx,
                                       const TrackingVolume& world,
                                       const Config& config,
                                       const Logger& logger) {
  if (!std::isfinite(config.fallbackWarningFraction) ||
      config.fallbackWarningFraction < 0 ||
      config.fallbackWarningFraction > 1 ||
      !std::isfinite(config.maxTolerance) || config.maxTolerance < 0 ||
      config.maxGridCells == 0) {
    throw std::invalid_argument("Invalid RZ volume lookup configuration");
  }
  m_impl = std::make_unique<Impl>(gctx, world, config);
  m_impl->warnAboutFallback(world, config.fallbackWarningFraction, logger);
}

RzGridVolumeLookup::~RzGridVolumeLookup() = default;

const RzGridVolumeLookup::Config& RzGridVolumeLookup::config() const {
  return m_impl->config;
}

const TrackingVolume* RzGridVolumeLookup::lowestTrackingVolume(
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

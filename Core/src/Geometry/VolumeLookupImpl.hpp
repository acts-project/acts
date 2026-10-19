// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Geometry/CylinderVolumeBounds.hpp"
#include "Acts/Geometry/TrackingVolume.hpp"
#include "Acts/Utilities/Logger.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <numbers>
#include <vector>

#include <Eigen/Geometry>

namespace Acts::detail {
// Select direct children linearly and descend through their configured policy.
// Accelerator fallbacks bypass dispatch on the root to avoid recursion cycles.
const TrackingVolume* linearVolumeScan(const TrackingVolume& root,
                                       const GeometryContext& gctx,
                                       const Vector3& position,
                                       double tolerance);

constexpr std::size_t noRzNode = std::numeric_limits<std::size_t>::max();

// Unlike Acts::AxisAlignedBoundingBox::intersect, Eigen's AlignedBox includes
// the upper edges, as required for conservative RZ candidate selection.
// Conservative envelope and hierarchy constraints for one subtree candidate.
struct RzLookupChild {
  const TrackingVolume* volume;
  Eigen::AlignedBox2d box;
  bool indexable = false;
  bool delegate = false;
  bool legacy = false;
  std::vector<const TrackingVolume*> ancestors;

  RzLookupChild(const GeometryContext& gctx, const TrackingVolume& v,
                double maxTolerance)
      : volume(&v) {
    const auto* cylinder =
        dynamic_cast<const CylinderVolumeBounds*>(&v.volumeBounds());
    const auto& trf = v.localToGlobalTransform(gctx);
    if (!v.isAlignable() && cylinder != nullptr && trf.linear().isIdentity(0) &&
        cylinder->get(CylinderVolumeBounds::eHalfPhiSector) ==
            std::numbers::pi) {
      const double offset = trf.translation().head<2>().norm();
      const double margin = maxTolerance + 1e-9;
      const double hz = cylinder->get(CylinderVolumeBounds::eHalfLengthZ);
      box.min() =
          Vector2(cylinder->get(CylinderVolumeBounds::eMinR) - offset - margin,
                  trf.translation().z() - hz - margin);
      box.max() =
          Vector2(cylinder->get(CylinderVolumeBounds::eMaxR) + offset + margin,
                  trf.translation().z() + hz + margin);
      indexable = true;
    }
  }
};

// One spatial index covers the subtree. Candidates are ordered in postorder:
// descendants before their parent fallback, and an entire earlier sibling's
// subtree before later siblings. Ancestor containment checks preserve the
// recursive scan even when child bounds extend outside their parents.
template <typename Partition>
struct RzLookupIndex {
  bool gen1 = false;
  bool gridLimit = false;
  bool hasIndex = false;
  std::size_t alignable = 0;
  std::size_t unsupported = 0;
  std::size_t legacy = 0;
  Eigen::AlignedBox2d domain;
  std::vector<RzLookupChild> children;
  std::vector<std::size_t> fallbackCandidates;
  Partition partition;

  RzLookupIndex(const GeometryContext& gctx, const TrackingVolume& root,
                const typename Partition::Config& cfg) {
    if (root.confinedVolumes() || !root.denseVolumes().empty()) {
      gen1 = true;
      return;
    }
    std::vector<const TrackingVolume*> ancestors;
    auto collect = [&](auto&& self, const TrackingVolume& volume) -> void {
      RzLookupChild entry(gctx, volume, cfg.maxTolerance);
      entry.ancestors = ancestors;
      entry.legacy =
          volume.volumeLookup() == nullptr &&
          (volume.confinedVolumes() || !volume.denseVolumes().empty());
      entry.delegate = volume.volumeLookup() != nullptr || entry.legacy;
      if (!entry.delegate) {
        ancestors.push_back(&volume);
        for (const auto& child : volume.volumes()) {
          self(self, child);
        }
        ancestors.pop_back();
      }
      if (entry.indexable) {
        domain.extend(entry.box);
        hasIndex = true;
      } else {
        fallbackCandidates.push_back(children.size());
        alignable += volume.isAlignable();
        unsupported += !volume.isAlignable();
      }
      legacy += entry.legacy;
      children.push_back(std::move(entry));
    };
    for (const auto& child : root.volumes()) {
      collect(collect, child);
    }
    if (hasIndex) {
      gridLimit = !partition.build(children, cfg);
    }
  }

  void warnAboutFallback(const TrackingVolume& volume, double fraction,
                         const Logger& constructionLogger) const {
    const auto logger = [&]() -> const Logger& { return constructionLogger; };
    if (fraction == 1 || !logger().doPrint(Logging::WARNING)) {
      return;
    }
    const auto total = gen1 ? std::size_t{1} : children.size();
    std::size_t affected = fallbackCandidates.size();
    // Indexed Gen1 boundaries still delegate to the legacy scan.
    for (const auto& child : children) {
      affected += child.indexable && child.legacy;
    }
    if (gen1 || gridLimit) {
      affected = total;
    }
    if (total == 0 || static_cast<double>(affected) / total <= fraction) {
      return;
    }
    ACTS_WARNING("RZ volume lookup for '"
                 << volume.volumeName() << "': " << affected << "/" << total
                 << " volumes in the indexed scope require linear fallback ("
                 << 100. * affected / total << "%, warning threshold "
                 << 100. * fraction << "%). Reasons: alignable=" << alignable
                 << ", unsupported shape/orientation=" << unsupported
                 << ", Gen1=" << (legacy + gen1)
                 << ", grid cell limit=" << gridLimit
                 << ". Explicit descendant lookup boundaries remain active.");
  }

  const TrackingVolume* find(const TrackingVolume& root,
                             const GeometryContext& gctx, const Vector3& p,
                             const Vector2& key, double tolerance) const {
    if (gen1 || gridLimit) {
      return linearVolumeScan(root, gctx, p, tolerance);
    }
    if (!root.inside(gctx, p, tolerance)) {
      return nullptr;
    }
    // IDs encode hierarchy priority. Ordered partitions can stop at the first
    // match; the existing point KDTree visits all matches in tree order.
    std::size_t best = noRzNode;
    auto consider = [&](std::size_t id) {
      if (id >= best) {
        return false;
      }
      const auto& entry = children[id];
      if ((entry.indexable && !entry.box.contains(key)) ||
          !entry.volume->inside(gctx, p, tolerance) ||
          !std::ranges::all_of(entry.ancestors,
                               [&](const TrackingVolume* ancestor) {
                                 return ancestor->inside(gctx, p, tolerance);
                               })) {
        return true;
      }
      best = id;
      return false;
    };
    for (auto id : fallbackCandidates) {
      if (!consider(id)) {
        break;
      }
    }
    if (hasIndex && domain.contains(key)) {
      partition.visitCandidates(key, consider);
    }
    if (best == noRzNode) {
      return &root;
    }
    const auto& entry = children[best];
    return entry.delegate
               ? entry.volume->lowestTrackingVolume(gctx, p, tolerance)
               : entry.volume;
  }
};
}  // namespace Acts::detail

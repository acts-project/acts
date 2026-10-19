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
#include "Acts/Utilities/VectorHelpers.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <numbers>
#include <stdexcept>
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
  // Stop indexing here and dispatch to this volume's explicitly set lookup.
  bool delegate = false;
  // Parents between this candidate and the lookup root, excluding both ends.
  // Flattening the hierarchy must not bypass their exact containment checks.
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
      const double offset = VectorHelpers::perp(trf.translation());
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
  bool gridLimit = false;
  bool hasIndex = false;
  std::size_t alignable = 0;
  std::size_t unsupported = 0;
  Eigen::AlignedBox2d domain;
  std::vector<RzLookupChild> children;
  std::vector<std::size_t> fallbackCandidates;
  Partition partition;

  RzLookupIndex(const GeometryContext& gctx, const TrackingVolume& root,
                const typename Partition::Config& cfg) {
    // A bare TrackingVolume has no generation flag. Reject Gen1-only storage
    // wherever it is encountered in this index's scope, including the root.
    auto requireGen3 = [](const TrackingVolume& volume) {
      if (volume.confinedVolumes() || volume.confinedLayers() != nullptr ||
          !volume.denseVolumes().empty()) {
        throw std::invalid_argument(
            "RZ volume lookup does not support Gen1 storage in volume '" +
            volume.volumeName() + "'");
      }
    };
    requireGen3(root);
    std::vector<const TrackingVolume*> ancestors;
    // Collect geometry independently of the spatial partition. Postorder IDs
    // encode recursive-scan priority; build() below indexes their RZ envelopes.
    auto collect = [&](auto&& self, const TrackingVolume& volume) -> void {
      requireGen3(volume);
      RzLookupChild entry(gctx, volume, cfg.maxTolerance);
      entry.ancestors = ancestors;
      entry.delegate = volume.volumeLookup() != nullptr;
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
        if (volume.isAlignable()) {
          ++alignable;
        } else {
          ++unsupported;
        }
      }
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
    const auto total = children.size();
    const auto affected = gridLimit ? total : fallbackCandidates.size();
    if (total == 0 ||
        static_cast<double>(affected) / static_cast<double>(total) <=
            fraction) {
      return;
    }
    ACTS_WARNING("RZ volume lookup for '"
                 << volume.volumeName() << "': " << affected << "/" << total
                 << " volumes in the indexed scope require linear fallback ("
                 << 100. * static_cast<double>(affected) /
                        static_cast<double>(total)
                 << "%, warning threshold " << 100. * fraction
                 << "%). Reasons: alignable=" << alignable
                 << ", unsupported shape/orientation=" << unsupported
                 << ", grid cell limit=" << gridLimit
                 << ". Explicit descendant lookup boundaries remain active.");
  }

  const TrackingVolume* find(const TrackingVolume& root,
                             const GeometryContext& gctx, const Vector3& p,
                             const Vector2& key, double tolerance) const {
    if (gridLimit) {
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

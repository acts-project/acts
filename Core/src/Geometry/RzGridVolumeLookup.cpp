// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Geometry/RzGridVolumeLookup.hpp"

#include "Acts/Geometry/CylinderVolumeBounds.hpp"
#include "Acts/Geometry/TrackingVolume.hpp"
#include "Acts/Utilities/Axis.hpp"
#include "Acts/Utilities/MultiAxis.hpp"
#include "Acts/Utilities/VectorHelpers.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <numbers>
#include <numeric>
#include <optional>
#include <span>
#include <stdexcept>
#include <vector>

#include <Eigen/Geometry>

#include "VolumeLookupImpl.hpp"

namespace Acts {
namespace {
constexpr std::size_t noRzNode = std::numeric_limits<std::size_t>::max();

// Unlike Acts::AxisAlignedBoundingBox::intersect, Eigen's AlignedBox includes
// the upper edges, as required for conservative RZ candidate selection.
// Conservative envelope for one subtree candidate.
struct RzVolumeCandidate {
  const TrackingVolume* m_volume;
  Eigen::AlignedBox2d m_box;
  bool m_indexable = false;
  // Stop indexing here and dispatch to this volume's explicitly set lookup.
  bool m_delegate = false;

  RzVolumeCandidate(const GeometryContext& gctx, const TrackingVolume& v,
                    double maxTolerance)
      : m_volume(&v) {
    const auto* cylinder =
        dynamic_cast<const CylinderVolumeBounds*>(&v.volumeBounds());
    const auto& trf = v.localToGlobalTransform(gctx);
    // Only fixed full cylinders with near-identity rotations are indexed.
    if (!v.isAlignable() && cylinder != nullptr &&
        trf.linear().isIdentity(1e-12) &&
        cylinder->get(CylinderVolumeBounds::eHalfPhiSector) ==
            std::numbers::pi) {
      const double offset = VectorHelpers::perp(trf.translation());
      const double hz = cylinder->get(CylinderVolumeBounds::eHalfLengthZ);
      // Bound displacement from the accepted rotation, including tolerance.
      // This keeps the envelope conservative even for tiny nonzero tilts.
      const double extent =
          std::hypot(cylinder->get(CylinderVolumeBounds::eMaxR) + maxTolerance,
                     hz + maxTolerance);
      const double margin =
          maxTolerance + 1e-9 +
          (trf.linear() - RotationMatrix3::Identity()).norm() * extent;
      m_box.min() =
          Vector2(cylinder->get(CylinderVolumeBounds::eMinR) - offset - margin,
                  trf.translation().z() - hz - margin);
      m_box.max() =
          Vector2(cylinder->get(CylinderVolumeBounds::eMaxR) + offset + margin,
                  trf.translation().z() + hz + margin);
      m_indexable = true;
    }
  }
};

// Grid edges are the sorted unique rectangle edges. Each rectangle occupies
// every intersected cell, including its upper-edge bin. Candidate lists retain
// hierarchy order and use CSR storage: offsets delimit slices in one array.
struct RzGrid {
  using Config = RzGridVolumeLookup::Config;
  using RzAxis = Axis<AxisType::Variable, AxisBoundaryType::Bound>;
  using RzAxes = MultiAxis<RzAxis, RzAxis>;
  std::optional<RzAxes> m_axes;
  std::vector<std::size_t> m_offsets;
  std::vector<std::size_t> m_candidates;

  bool build(const std::vector<RzVolumeCandidate>& children,
             const Config& config) {
    std::array<std::vector<double>, 2> edges;
    for (int axis = 0; axis < 2; ++axis) {
      auto& e = edges[axis];
      for (const auto& child : children) {
        if (!child.m_indexable) {
          continue;
        }
        e.push_back(child.m_box.min()[axis]);
        e.push_back(child.m_box.max()[axis]);
      }
      std::ranges::sort(e);
      const auto duplicates = std::ranges::unique(e);
      e.erase(duplicates.begin(), duplicates.end());
    }
    // N edges give N-1 regular bins plus two reserved under-/overflow bins.
    const auto nr = edges[0].size() + 1;
    const auto nz = edges[1].size() + 1;
    // Check nr*nz against the limit without overflowing the product.
    if (nr > config.maxGridCells / nz) {
      return false;
    }
    m_axes.emplace(RzAxis(std::move(edges[0])), RzAxis(std::move(edges[1])));
    // One extra offset marks the end of the last cell's candidate slice.
    m_offsets.resize(m_axes->getNTotalBins(true) + 1);
    auto cells = [&](const Eigen::AlignedBox2d& b, auto&& visit) {
      const auto lower = m_axes->getLocalBinsFromPoint(b.min());
      const auto upper = m_axes->getLocalBinsFromPoint(b.max());
      // Include the bin at the upper edge. Bound axes clamp the final edge to
      // the last regular bin, matching the inclusive rectangle containment.
      for (auto r = lower[0]; r <= upper[0]; ++r) {
        for (auto z = lower[1]; z <= upper[1]; ++z) {
          visit(m_axes->getGlobalBinFromLocalBins({r, z}));
        }
      }
    };
    // Shift counts by one, prefix-sum to slice offsets, then fill with separate
    // cursors so the offsets remain unchanged.
    for (const auto& child : children) {
      if (!child.m_indexable) {
        continue;
      }
      cells(child.m_box, [&](auto cell) { ++m_offsets[cell + 1]; });
    }
    std::partial_sum(m_offsets.begin(), m_offsets.end(), m_offsets.begin());
    m_candidates.resize(m_offsets.back());
    auto cursor = m_offsets;
    for (std::size_t i = 0; i < children.size(); ++i) {
      if (!children[i].m_indexable) {
        continue;
      }
      cells(children[i].m_box,
            [&](auto cell) { m_candidates[cursor[cell]++] = i; });
    }
    return true;
  }

  std::span<const std::size_t> candidates(const Vector2& key) const {
    // The caller checks the inclusive domain; MultiAxis::isInside() would
    // exclude the upper edge, while Bound axes alone also clamp outside points.
    const auto cell = m_axes->getGlobalBinFromPoint(key);
    return std::span{m_candidates}.subspan(
        m_offsets[cell], m_offsets[cell + 1] - m_offsets[cell]);
  }
};
}  // namespace
// One spatial index covers the subtree. Candidates are ordered in postorder:
// descendants before their parent fallback, and an entire earlier sibling's
// subtree before later siblings. Geometry construction must ensure each child
// is contained within its parent, so candidates need no ancestor checks.
struct RzGridVolumeLookup::Impl {
  const TrackingVolume& m_root;
  Config m_config;
  bool m_gridLimit = false;
  bool m_hasIndex = false;
  std::size_t m_alignable = 0;
  std::size_t m_unsupported = 0;
  Eigen::AlignedBox2d m_domain;
  std::vector<RzVolumeCandidate> m_children;
  std::vector<std::size_t> m_fallbackCandidates;
  // One grid for this lookup's entire collected subtree, not one per volume.
  RzGrid m_grid;

  Impl(const GeometryContext& gctx, const TrackingVolume& root,
       const Config& cfg)
      : m_root(root), m_config(cfg) {
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
    // Collect geometry independently of the spatial grid. Postorder IDs
    // encode recursive-scan priority; build() below indexes their RZ envelopes.
    auto collect = [&](auto&& self, const TrackingVolume& volume) -> void {
      requireGen3(volume);
      RzVolumeCandidate entry(gctx, volume, cfg.maxTolerance);
      entry.m_delegate = volume.volumeLookup() != nullptr;
      if (!entry.m_delegate) {
        for (const auto& child : volume.volumes()) {
          self(self, child);
        }
      }
      if (entry.m_indexable) {
        m_domain.extend(entry.m_box);
        m_hasIndex = true;
      } else {
        m_fallbackCandidates.push_back(m_children.size());
        if (volume.isAlignable()) {
          ++m_alignable;
        } else {
          ++m_unsupported;
        }
      }
      m_children.push_back(std::move(entry));
    };
    for (const auto& child : root.volumes()) {
      collect(collect, child);
    }
    // Build the single grid after collecting the whole indexed scope.
    if (m_hasIndex) {
      m_gridLimit = !m_grid.build(m_children, cfg);
    }
  }

  void warnAboutFallback(const TrackingVolume& volume, double fraction,
                         const Logger& logger) const {
    // A threshold of 1 suppresses this warning.
    if (fraction == 1 || !logger.doPrint(Logging::WARNING)) {
      return;
    }
    const auto total = m_children.size();
    const auto affected = m_gridLimit ? total : m_fallbackCandidates.size();
    // Avoid division by zero; warn only above the configured fraction.
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
                 << "%). Reasons: alignable=" << m_alignable
                 << ", unsupported shape/orientation=" << m_unsupported
                 << ", grid cell limit=" << m_gridLimit
                 << ". Explicit descendant lookup boundaries remain active.");
  }

  const TrackingVolume* find(const GeometryContext& gctx, const Vector3& p,
                             const Vector2& key, double tolerance) const {
    if (m_gridLimit) {
      return detail::linearVolumeScan(m_root, gctx, p, tolerance);
    }
    if (!m_root.inside(gctx, p, tolerance)) {
      return nullptr;
    }
    // IDs encode hierarchy priority. Each ordered candidate list can stop
    // at its first exact match; choose the earlier match across both lists.
    std::size_t best = noRzNode;
    // Return true to keep scanning this ordered list, false to stop it.
    auto consider = [&](std::size_t id) {
      // Later IDs cannot beat a match already found in either list.
      if (id >= best) {
        return false;
      }
      const auto& entry = m_children[id];
      // Check the conservative box when available, then exact containment.
      // A miss continues the search; a match stops this ordered list.
      if ((entry.m_indexable && !entry.m_box.contains(key)) ||
          !entry.m_volume->inside(gctx, p, tolerance)) {
        return true;
      }
      best = id;
      // The other list may still contain a higher-priority match.
      return false;
    };
    for (auto id : m_fallbackCandidates) {
      if (!consider(id)) {
        break;
      }
    }
    if (m_hasIndex && m_domain.contains(key)) {
      for (auto id : m_grid.candidates(key)) {
        if (!consider(id)) {
          break;
        }
      }
    }
    if (best == noRzNode) {
      return &m_root;
    }
    const auto& entry = m_children[best];
    return entry.m_delegate
               ? entry.m_volume->lowestTrackingVolume(gctx, p, tolerance)
               : entry.m_volume;
  }
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
  return m_impl->m_config;
}

const TrackingVolume* RzGridVolumeLookup::findLowest(
    const GeometryContext& gctx, const Vector3& position,
    double tolerance) const {
  if (!std::isfinite(tolerance) || tolerance < 0 ||
      tolerance > m_impl->m_config.maxTolerance) {
    return detail::linearVolumeScan(m_impl->m_root, gctx, position, tolerance);
  }
  return m_impl->find(gctx, position,
                      Vector2(VectorHelpers::perp(position), position.z()),
                      tolerance);
}
}  // namespace Acts

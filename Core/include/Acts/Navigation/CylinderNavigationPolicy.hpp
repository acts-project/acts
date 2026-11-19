// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Navigation/INavigationPolicy.hpp"
#include "Acts/Navigation/NavigationStream.hpp"

namespace Acts {

class TrackingVolume;
class GeometryContext;
class Logger;

/// Optimized navigation policy for cylindrical volumes that intelligently
/// selects which portals to add as candidates based on geometric analysis.
///
/// This policy performs geometric calculations to determine which cylinder
/// faces (inner/outer cylinder, positive/negative discs) are reachable from
/// the current position and direction, avoiding unnecessary intersection
/// calculations with unreachable portals.
///
/// Algorithm overview:
/// 1. Transform position/direction to volume-local coordinates
/// 2. Check disc intersection: If not parallel to the disc, calculate
///    intersection with appropriate disc (positive/negative based on
///    direction) and verify if intersection point lies within the annular disc
///    bounds
/// 3. Check inner cylinder intersection: Use transverse dot and cross
///    products to test whether the forward ray passes inside the inner
///    radius, without normalizing the direction or solving a quadratic
/// 4. Optimization: If a disc is hit, require the closest approach to the
///    z-axis to precede it. An inner cylinder hit then blocks the disc,
///    which can be discarded
/// 5. Default fallback: If neither disc nor inner cylinder are hit, add outer
///    cylinder as the target (particle will exit through outer boundary)
///
/// Constraints:
/// - Only works with full-phi cylindrical volumes having non-zero inner radius
/// - Requires exactly 4 portals (inner cylinder, outer cylinder, ±z discs),
///   which also means the volume cannot confine sub-volumes, as it would carry
///   their portals as well
/// - Only produces portal candidates: a volume confining surfaces needs to
///   combine this policy with one that adds the surface candidates
///
/// Performance benefits:
/// - Reduces intersection calculations by ~2-3x compared to trying all portals
/// - Eliminates false positive intersections that would be filtered later
class CylinderNavigationPolicy final : public INavigationPolicy {
 public:
  /// Check if this policy is applicable to a volume, i.e. if the volume
  /// satisfies all of the constraints documented above. Reasons for a volume
  /// being rejected are reported to the logger.
  /// @param volume is the volume to check
  /// @param logger is the logger
  /// @return true if the policy can be constructed from @p volume
  static bool isApplicable(const TrackingVolume& volume, const Logger& logger);

  /// Constructor from a volume
  /// @note The volume needs to satisfy @ref isApplicable, otherwise the
  ///       constructor throws @c std::invalid_argument
  /// @param gctx is the geometry context
  /// @param volume is the volume to navigate
  /// @param logger is the logger
  CylinderNavigationPolicy(const GeometryContext& gctx,
                           const TrackingVolume& volume, const Logger& logger);

  /// Intelligently select and add portal candidates based on geometric analysis
  ///
  /// The algorithm determines which portals are geometrically reachable:
  /// - Disc portals: Added if ray intersects the annular disc area
  /// - Inner cylinder: Added if ray's closest approach to z-axis is within
  /// inner radius
  /// - Outer cylinder: Added as fallback when no other portals are reachable
  ///
  /// @param gctx is the geometry context
  /// @param args are the navigation arguments containing position and direction
  /// @param state is the navigation policy state
  /// @param stream is the navigation stream to update with selected candidates
  /// @param logger is the logger for debugging output
  void initializeCandidates(const GeometryContext& gctx,
                            const NavigationArguments& args,
                            NavigationPolicyState& state,
                            AppendOnlyNavigationStream& stream,
                            const Logger& logger) const;

  /// Connect the policy to a navigation delegate
  /// @param delegate is the navigation delegate
  void connect(NavigationDelegate& delegate) const override;

 private:
  /// Pointer to the cylindrical tracking volume this policy navigates
  const TrackingVolume* m_volume;

  /// Inverse transform for converting global to local coordinates (optional)
  /// Only computed if volume transform is not identity
  std::optional<Transform3> m_itransform;

  /// Cached cylinder bounds for efficient access during navigation
  double m_halfLengthZ;  ///< Half-length of cylinder in z-direction
  double m_rMin2;        ///< Squared inner radius for fast comparison
  double m_rMax2;        ///< Squared outer radius for fast comparison

  /// Direct portal references for efficient candidate addition
  /// Index corresponds to CylinderVolumeBounds::Face enum values
  std::array<const Portal*, 4> m_portals{};
};

static_assert(NavigationPolicyConcept<CylinderNavigationPolicy>);

}  // namespace Acts

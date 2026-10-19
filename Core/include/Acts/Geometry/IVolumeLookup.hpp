// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Geometry/GeometryContext.hpp"

#include <functional>
#include <memory>

namespace Acts {
class Logger;
class TrackingVolume;

/// Immutable point-to-volume lookup policy owned by a TrackingVolume.
/// Implementations must preserve lowestTrackingVolume containment, tolerance,
/// first-matching-sibling precedence and parent fallback semantics. Referenced
/// volumes must outlive the policy and retain their addresses. Queries must
/// support concurrent readers; install policies before concurrent use.
class IVolumeLookup {
 public:
  virtual ~IVolumeLookup() = default;

  /// Find the deepest containing volume, or nullptr outside the bound root
  /// volume. Respect explicit descendant lookups as delegation boundaries.
  /// @param gctx Geometry context for containment checks
  /// @param position Global query position
  /// @param tolerance Containment tolerance
  /// @return Containing tracking volume
  virtual const TrackingVolume* lowestTrackingVolume(
      const GeometryContext& gctx, const Vector3& position,
      double tolerance) const = 0;
};
/// Construct a non-null lookup for a volume after its descendants are
/// finalized.
using VolumeLookupFactory = std::function<std::unique_ptr<IVolumeLookup>(
    const GeometryContext&, const TrackingVolume&, const Logger&)>;
}  // namespace Acts

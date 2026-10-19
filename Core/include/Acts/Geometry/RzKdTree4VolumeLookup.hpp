// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Geometry/IVolumeLookup.hpp"
#include "Acts/Utilities/Logger.hpp"

#include <cstddef>
#include <memory>

namespace Acts {
/// RZ lookup using Acts::KDTree with (rmin, zmin, rmax, zmax) point keys.
/// Uses the existing tree's fixed leaf size of four.
/// Rebuild after changing
/// volume bounds, transforms, hierarchy or descendant lookup policies. Indexes
/// the full subtree, stopping at explicitly configured descendant lookups.
/// Unsupported bounds are checked linearly alongside indexed candidates.
class RzKdTree4VolumeLookup final : public IVolumeLookup {
 public:
  /// Serializable construction parameters; assignments are rebuilt from
  /// volumes.
  struct Config {
    /// Largest accelerated tolerance; other tolerances use linear lookup
    double maxTolerance = 1e-3;
    /// Warn when more than this fraction of volumes in the indexed scope
    /// requires linear fallback. Explicit descendant policies remain active.
    /// Range [0, 1]; 1 disables the warning. Not a query-frequency metric.
    double fallbackWarningFraction = 0.1;
  };

  /// Build an immutable index over a finalized subtree.
  /// @param gctx Construction context
  /// @param world Root of the indexed subtree
  /// @param config Construction settings
  /// @param logger Construction diagnostics
  RzKdTree4VolumeLookup(const GeometryContext& gctx,
                        const TrackingVolume& world, const Config& config,
                        const Logger& logger = getDummyLogger());
  ~RzKdTree4VolumeLookup() override;

  /// @copydoc IVolumeLookup::lowestTrackingVolume
  const TrackingVolume* lowestTrackingVolume(const TrackingVolume& root,
                                             const GeometryContext& gctx,
                                             const Vector3& position,
                                             double tolerance) const override;

  /// Return construction parameters for serialization.
  /// @return Construction settings
  const Config& config() const;

 private:
  struct Impl;
  std::unique_ptr<const Impl> m_impl;
};
}  // namespace Acts

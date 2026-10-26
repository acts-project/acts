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
/// Built-in RZ grid spatial index. Construct after geometry closure, once
/// bounds, hierarchy and descendant policies are finalized; prefer installing
/// through a Blueprint or StaticBlueprintNode lookup factory.
/// Rebuild after changing
/// volume bounds, transforms, hierarchy or descendant lookup policies. Indexes
/// the full subtree, stopping at explicitly configured descendant lookups.
/// Child volumes must be contained within their parents.
/// Unsupported bounds are checked linearly alongside indexed candidates.
class RzGridVolumeLookup final : public IVolumeLookup {
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
    /// Limit allocated grid cells, including reserved under-/overflow slots;
    /// larger grids use linear lookup
    std::size_t maxGridCells = 1000000;
  };

  /// Build an immutable index over a finalized subtree.
  /// @param gctx Construction context
  /// @param world Root of the indexed subtree
  /// @param config Construction settings
  /// @param logger Construction diagnostics
  /// @throws std::invalid_argument For invalid configuration or Gen1 storage
  /// (confined layers/volumes or dense volumes) encountered in the indexed
  /// scope.
  RzGridVolumeLookup(const GeometryContext& gctx, const TrackingVolume& world,
                     const Config& config,
                     const Logger& logger = getDummyLogger());
  ~RzGridVolumeLookup() override;

  /// @copydoc IVolumeLookup::findLowest
  const TrackingVolume* findLowest(const GeometryContext& gctx,
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

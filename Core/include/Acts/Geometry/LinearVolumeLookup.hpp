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

namespace Acts {
/// Default lookup policy using the volume hierarchy's existing recursive scan.
class LinearVolumeLookup final : public IVolumeLookup {
 public:
  /// Construction settings.
  struct Config {
    /// Warn above this sibling count within the scanned subtree; zero disables
    /// warnings.
    std::size_t warningThreshold = 64;
  };

  /// Construct a lookup and diagnose a finalized subtree.
  /// @param world Root of the scanned subtree
  /// @param config Construction settings
  /// @param logger Construction diagnostics
  LinearVolumeLookup(const TrackingVolume& world, const Config& config,
                     const Logger& logger = getDummyLogger());

  /// Return construction settings, including the warning threshold.
  /// @return Construction settings
  const Config& config() const;

  /// @copydoc IVolumeLookup::lowestTrackingVolume
  const TrackingVolume* lowestTrackingVolume(const GeometryContext& gctx,
                                             const Vector3& position,
                                             double tolerance) const override;

 private:
  void warnIfNeeded(const Logger& logger) const;

  const TrackingVolume& m_root;
  Config m_config;
};
}  // namespace Acts

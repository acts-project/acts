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

/// Policy which adds **all** candidates of the configured type to the
/// stream
class TryAllNavigationPolicy final : public INavigationPolicy {
 public:
  struct Config {
    bool portals = true;
    bool sensitives = true;
  };

  /// Constructor from a volume
  /// @param config The configuration for the policy
  /// @param gctx is the geometry context
  /// @param volume is the volume to navigate
  /// @param logger is the logger
  TryAllNavigationPolicy(const Config& config, const GeometryContext& gctx,
                         const TrackingVolume& volume, const Logger& logger);

  TryAllNavigationPolicy(const GeometryContext& gctx,
                         const TrackingVolume& volume, const Logger& logger);

  /// Update the navigation state with all volume portals.
  /// @param args are the navigation arguments
  void initializeCandidates(const NavigationArguments& args,
                            AppendOnlyNavigationStream& stream,
                            const Logger& logger) const;

  /// Connect the policy to a navigation delegate
  /// @param delegate is the navigation delegate
  void connect(NavigationDelegate& delegate) const override;

 private:
  Config m_cfg;
  const TrackingVolume* m_volume;
};

static_assert(NavigationPolicyConcept<TryAllNavigationPolicy>);

}  // namespace Acts

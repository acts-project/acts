// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Navigation/INavigationPolicy.hpp"

namespace Acts {

class TrackingVolume;
class GeometryContext;
class Logger;

/// Policy which adds **all** portals of the volume to the main navigation
/// stream
class TryAllPortalNavigationPolicy final : public INavigationPolicy {
 public:
  /// Constructor from a volume
  /// @param gctx is the geometry context
  /// @param volume is the volume to navigate
  /// @param logger is the logger
  explicit TryAllPortalNavigationPolicy(const GeometryContext& gctx,
                                        const TrackingVolume& volume,
                                        const Logger& logger);

  /// Update the navigation state with all volume portals.
  /// @param args are the navigation arguments
  void updateState(const NavigationArguments& args) const;

  /// Connect the policy to a navigation delegate
  /// @param delegate is the navigation delegate
  void connect(NavigationDelegate& delegate) const override;

 private:
  const TrackingVolume* m_volume;
};

static_assert(NavigationPolicyConcept<TryAllPortalNavigationPolicy>);

/// Policy which adds **all** surfaces of the volume to the main navigation
/// @note This does **not** include portal surfaces, but does not make any distinction
///       between sensitive and other surfaces in the volume.
class TryAllSurfaceNavigationPolicy final : public INavigationPolicy {
 public:
  /// Constructor from a volume
  /// @param gctx is the geometry context
  /// @param volume is the volume to navigate
  /// @param logger is the logger
  explicit TryAllSurfaceNavigationPolicy(const GeometryContext& gctx,
                                         const TrackingVolume& volume,
                                         const Logger& logger);

  /// Update the navigation state with all volume surfaces.
  /// @param args are the navigation arguments
  void updateState(const NavigationArguments& args) const;

  /// Connect the policy to a navigation delegate
  /// @param delegate is the navigation delegate
  void connect(NavigationDelegate& delegate) const override;

 private:
  const TrackingVolume* m_volume;
};

static_assert(NavigationPolicyConcept<TryAllSurfaceNavigationPolicy>);

}  // namespace Acts

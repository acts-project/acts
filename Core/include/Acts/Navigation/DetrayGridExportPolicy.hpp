// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Units.hpp"
#include "Acts/Navigation/INavigationPolicy.hpp"

#include <span>

namespace Acts {

class TrackingVolume;
class GeometryContext;
class Logger;
class SurfaceArrayNavigationPolicy;

/// Navigation policy that describes how the surface grid of its volume is
/// exported to detray. It does not take part in ACTS navigation.
///
/// It refers to the @ref SurfaceArrayNavigationPolicy constructed before it
/// by the same @ref NavigationPolicyFactory, and carries the track family
/// the exported grid must serve. The Detray plugin uses it to fill each bin of
/// the detray grid with the surfaces the ACTS lookup returns for these tracks,
/// so detray can navigate the grid without a search window.
class DetrayGridExportPolicy final : public INavigationPolicy {
 public:
  /// Track family the exported grid must serve
  struct Config {
    /// Half range of the vertex z positions around z = 0
    double vertexZHalfRange = 200. * UnitConstants::mm;
    /// Lowest transverse momentum whose bending is covered. A value <= 0
    /// only covers straight lines.
    double minPt = 1. * UnitConstants::GeV;
    /// Solenoid field strength used to estimate the bending
    double bField = 2. * UnitConstants::T;
    /// Fill every bin with its maximum ACTS neighbor window instead of the
    /// track-dependent window. Independent of the track family, so it also
    /// serves displaced tracks, at the cost of larger bins.
    bool fixedWindow = false;
  };

  /// Constructor
  /// @param gctx The geometry context
  /// @param volume The volume this policy belongs to
  /// @param logger The logger
  /// @param siblings Policies constructed before this one by the same
  ///        factory; exactly one @ref SurfaceArrayNavigationPolicy must be
  ///        among them, possibly nested in a composite policy
  /// @param config The configuration
  /// @throws std::invalid_argument if no or more than one
  ///         @ref SurfaceArrayNavigationPolicy is found
  DetrayGridExportPolicy(const GeometryContext& gctx,
                         const TrackingVolume& volume, const Logger& logger,
                         std::span<const INavigationPolicy* const> siblings,
                         const Config& config);

  /// Does nothing: the policy only carries export information
  /// @param gctx The geometry context
  /// @param args The navigation arguments
  /// @param state The navigation policy state
  /// @param stream The navigation stream, left unchanged
  /// @param logger The logger
  void initializeCandidates(const GeometryContext& gctx,
                            const NavigationArguments& args,
                            NavigationPolicyState& state,
                            AppendOnlyNavigationStream& stream,
                            const Logger& logger) const;

  /// Connect the no-op candidate initialization to a delegate
  /// @param delegate The navigation delegate to connect to
  void connect(NavigationDelegate& delegate) const override;

  /// Constant access to the configuration
  /// @return The configuration
  const Config& config() const;

  /// The sibling policy whose surface array is exported
  /// @return The surface array navigation policy
  const SurfaceArrayNavigationPolicy& surfaceArrayPolicy() const;

 private:
  Config m_cfg;
  // Owned by the same MultiNavigationPolicy, so it outlives this policy's use
  const SurfaceArrayNavigationPolicy* m_surfaceArrayPolicy = nullptr;
};

static_assert(NavigationPolicyConcept<DetrayGridExportPolicy>);

}  // namespace Acts

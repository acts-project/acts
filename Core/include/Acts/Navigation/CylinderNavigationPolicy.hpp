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

/// Policy which adds only the outer portals of a cylinder. It uses some
/// simplified calculations which set of portals to add to the candidate vector.
/// This policy does not handle contained volumes at all.
class CylinderNavigationPolicy final : public INavigationPolicy {
 public:
  struct Config {};

  /// Constructor from a volume
  /// @param gctx is the geometry context
  /// @param volume is the volume to navigate
  /// @param logger is the logger
  /// @param config The configuration for the policy
  CylinderNavigationPolicy(const GeometryContext& gctx,
                           const TrackingVolume& volume, const Logger& logger,
                           const Config& config);

  /// Add all candidates to the stream
  /// @param args are the navigation arguments
  /// @param stream is the navigation stream to update
  /// @param logger is the logger
  void initializeCandidates(const NavigationArguments& args,
                            AppendOnlyNavigationStream& stream,
                            const Logger& logger) const;

  /// Connect the policy to a navigation delegate
  /// @param delegate is the navigation delegate
  void connect(NavigationDelegate& delegate) const override;

 private:
  Config m_cfg;
  const TrackingVolume* m_volume;
  std::optional<Transform3> m_itransform;

  double m_halfLengthZ;
  double m_rMin2;
  double m_rMax2;

  std::array<const Portal*, 4> m_portals;
};

static_assert(NavigationPolicyConcept<CylinderNavigationPolicy>);

}  // namespace Acts

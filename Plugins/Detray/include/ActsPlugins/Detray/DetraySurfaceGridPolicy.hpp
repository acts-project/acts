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
#include "Acts/Navigation/NavigationStream.hpp"

#include <cstddef>
#include <optional>

namespace Acts {
class GeometryContext;
class Logger;
class TrackingVolume;
}  // namespace Acts

namespace ActsPlugins {

/// Parameters of the detray surface grid of a volume.
///
/// Detray looks up a surface grid by intersecting the tangent of the track
/// with a virtual reference surface and reading a single bin. The grid is
/// filled so that this bin lists every surface a track of the configured
/// family can reach from there.
struct DetraySurfaceGridConfig {
  /// Half range of the vertex z positions around z = 0
  double vertexZHalfRange = 150. * Acts::UnitConstants::mm;
  /// Lowest transverse momentum covered; <= 0 covers straight tracks only
  double minPt = 0.5 * Acts::UnitConstants::GeV;
  /// Solenoid field strength used for the bending
  double bField = 2. * Acts::UnitConstants::T;
  /// Number of vertex z values sampled across the range
  std::size_t nVertexZ = 5;
  /// Number of curvature values sampled across [-max, max]
  std::size_t nCurvature = 5;
  /// Maximum spacing of the outline samples, in bins
  double maxSampleSpacing = 0.25;
  /// Bin count factor for axes the surface array bins automatically, relative
  /// to one bin per module. If unset, the surface array's own binning is kept.
  std::optional<double> numberOfBinsFactor = 2.;
};

/// Navigation policy carrying the @ref DetraySurfaceGridConfig of a volume.
///
/// It does not take part in ACTS navigation. When the volume also has a
/// surface array, the detray conversion builds the volume's detray surface
/// grid with these parameters instead of the converter defaults.
class DetraySurfaceGridPolicy final : public Acts::INavigationPolicy {
 public:
  /// Constructor
  /// @param gctx The geometry context
  /// @param volume The volume this policy belongs to
  /// @param logger The logger
  /// @param config The grid parameters
  DetraySurfaceGridPolicy(const Acts::GeometryContext& gctx,
                          const Acts::TrackingVolume& volume,
                          const Acts::Logger& logger,
                          const DetraySurfaceGridConfig& config);

  /// Constructor with default parameters
  /// @param gctx The geometry context
  /// @param volume The volume this policy belongs to
  /// @param logger The logger
  DetraySurfaceGridPolicy(const Acts::GeometryContext& gctx,
                          const Acts::TrackingVolume& volume,
                          const Acts::Logger& logger);

  /// Adds no candidates
  /// @param gctx The geometry context
  /// @param args The navigation arguments
  /// @param state The navigation policy state
  /// @param stream The navigation stream, left unchanged
  /// @param logger The logger
  void initializeCandidates(const Acts::GeometryContext& gctx,
                            const Acts::NavigationArguments& args,
                            Acts::NavigationPolicyState& state,
                            Acts::AppendOnlyNavigationStream& stream,
                            const Acts::Logger& logger) const;

  /// Connect the no-op candidate initialization to a delegate
  /// @param delegate The navigation delegate
  void connect(Acts::NavigationDelegate& delegate) const override;

  /// The grid parameters
  /// @return The configuration
  const DetraySurfaceGridConfig& config() const;

 private:
  DetraySurfaceGridConfig m_cfg;
};

static_assert(Acts::NavigationPolicyConcept<DetraySurfaceGridPolicy>);

}  // namespace ActsPlugins

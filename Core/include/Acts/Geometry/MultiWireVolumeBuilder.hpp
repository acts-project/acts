// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Geometry/NavigationPolicyFactory.hpp"
#include "Acts/Geometry/TrackingVolume.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/AxisSpec.hpp"
#include "Acts/Utilities/Logger.hpp"

namespace Acts {

/// @class MultiWireVolumeBuilder
/// @brief A class to build multiwire tracking volumes (e.g wire chambers)
class MultiWireVolumeBuilder {
 public:
  /// The axis configuration for the binning: axis direction, plus the bin
  /// expansion
  using Binning = std::tuple<AxisDirection, std::size_t>;

  /// Configuration Struct
  struct Config {
    /// The name of the tracking volume
    std::string name = "undefined";

    /// The surfaces to be wrapped from the tracking volume
    std::vector<std::shared_ptr<Surface>> mlSurfaces{};

    /// The local -> global transform of the tracking volume
    Transform3 transform{Transform3::Identity()};

    /// Connect the tracking geometry with an alignable volume placement
    /// Used instead of the transform if set
    VolumePlacementBase* alignablePlacement{};

    /// The bounds of the tracking volume
    std::shared_ptr<Acts::VolumeBounds> bounds = nullptr;

    /// Binning configuration for multi-wire volume
    std::vector<Binning> binning{};

    /// Boolean flag if staggering corrections should be applied
    bool correctOffsets{true};

    // The direction of the axis the shift of the surfaces is applied
    // It is the direction along the tubes
    AxisDirection shiftDirection{};
  };

  /// Constructor
  /// @param config The configuration struct
  /// @param logger The logger instance for screen output
  explicit MultiWireVolumeBuilder(
      const Config& config,
      std::unique_ptr<const Acts::Logger> logger = Acts::getDefaultLogger(
          "MultiWireVolumeBuilder", Acts::Logging::INFO));

  /// @brief Constructs the tracking volume with the wrapped surfaces
  /// @return a unique ptr of the tracking volume
  std::unique_ptr<Acts::TrackingVolume> buildVolume() const;

  /// @brief Creates a multilayer navigation policy factory that can be used for the trackingVolume
  /// or attached to a blueprint node
  /// @return Unique pointer to the created navigation policy factory
  std::unique_ptr<NavigationPolicyFactory> createNavigationPolicyFactory(
      const GeometryContext& gctx) const;

 private:
  // The config
  Config m_config;

  // The ACTS logger
  const Acts::Logger& logger() const { return *m_logger; }

  std::unique_ptr<const Acts::Logger> m_logger;

  /// Helper method to derive the grid parameters (both axes and the per-layer
  /// staggering correction) from the tube surfaces.
  ///
  /// The tube centers are projected, in the volume-local frame, onto the two
  /// binning directions. The shift axis (m_config.shiftDirection) is the one
  /// the tubes are staggered along; the layer axis is the other. Bins are
  /// centered on the tube lattice (edges at the inter-tube midpoints). The
  /// returned shift vector is indexed by layer-axis bin and corrects the
  /// per-layer stagger during bin registration.
  ///
  /// @param gctx the geometry context
  /// @return {shiftAxis, layerAxis, layerShifts} — layerShifts has one entry
  ///         per layer-axis bin.
  std::tuple<AxisSpec, AxisSpec, std::vector<double>> deriveGridParameters(
      const GeometryContext& gctx) const;
};

namespace Experimental {
/// @deprecated The blueprint geometry moved out of the `Acts::Experimental`
///             namespace. Use @ref Acts::MultiWireVolumeBuilder instead. This
///             alias is kept for backward compatibility and will be removed.
using MultiWireVolumeBuilder
    [[deprecated("Acts::Experimental::MultiWireVolumeBuilder moved to "
                 "Acts::MultiWireVolumeBuilder")]] =
        Acts::MultiWireVolumeBuilder;
}  // namespace Experimental

}  // namespace Acts

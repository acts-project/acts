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
#include "Acts/Utilities/Logger.hpp"
#include "Acts/Utilities/ProtoAxis.hpp"

namespace Acts::Experimental {

/// @class MultiWireVolumeBuilder
/// @brief A class to build multiwire tracking volumes (e.g wire chambers)
class MultiWireVolumeBuilder {
 public:
  /// The axis configuration for the binning
  using Binning = std::tuple<DirectedProtoAxis, std::size_t>;
  /// Configuration Struct
  struct Config {
    /// The name of the tracking volume
    std::string name = "undefined";

    /// The surfaces to be wrapped from the tracking volume
    std::vector<std::shared_ptr<Surface>> mlSurfaces = {};

    /// The transform of the tracking volume
    Transform3 transform = Transform3::Identity();

    /// The bounds of the tracking volume
    std::shared_ptr<Acts::VolumeBounds> bounds = nullptr;

    /// Binning configuration for multi-wire volume
    std::vector<Binning> binning = {};
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
  std::unique_ptr<NavigationPolicyFactory> createNavigationPolicyFactory()
      const;

 private:
  Config m_config;

  const Acts::Logger& logger() const { return *m_logger; }

  std::unique_ptr<const Acts::Logger> m_logger;
};

}  // namespace Acts::Experimental

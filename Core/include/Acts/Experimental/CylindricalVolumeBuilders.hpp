// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Geometry/Extent.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Utilities/Logger.hpp"

#include <exception>
#include <memory>
#include <string>
#include <vector>

/// @note This file is foreseen for the `Geometry` module

namespace Acts {

class DetectorVolume;
class InternalBlueprint;

struct CylindricalVolumeBuilder {
  /// Constructor with default logger
  ///
  /// @param level is the logging level
  CylindricalVolumeBuilder(Logging::Level level = Logging::INFO)
      : m_logger(getDefaultLogger("CylindricalVolumeBuilder", level)) {}

  /// Call operator that creates a single cylindrical volume
  ///
  /// @param gctx the geometry context for this building call
  /// @param iBlueprints the blue prints of the internal volume structure
  /// @param restriction is an extend that restricts the volume
  /// @param name the base name of the volumes
  ///
  /// @note throws an exception if more than one Internal blueprints
  /// are provided - and if the internal blueprint does not fit into
  /// the external restriction
  ///
  /// @return the new (cylindrical) container volume
  std::vector<std::shared_ptr<DetectorVolume>> operator()(
      const GeometryContext& gctx,
      const std::vector<InternalBlueprint>& iBlueprints,
      const Extent& restriction = Extent(false),
      const std::string& name = "Volume") noexcept(false);

  /// Private access method to the logger
  const Logger& logger() const { return *m_logger; }

  /// the logging instance
  std::unique_ptr<const Logger> m_logger;
};

struct CylindricalVolumeBuilderR {
  /// Constructor with default logger
  ///
  /// @param level is the logging level
  CylindricalVolumeBuilderR(Logging::Level level = Logging::INFO)
      : m_logger(getDefaultLogger("CylindricalVolumeBuilderR", level)) {}

  /// Call operator that creates volumes for a
  /// cylindrical container ordered in R
  ///
  /// @param gctx the geometry context for this building call
  /// @param iBlueprints the blue prints of the internal volume structure
  /// @param restriction is an extend that restricts the volume
  /// @param name the base name of the volumes
  ///
  /// In case no internal or only one internal blueprint is handed over
  /// this function will attempt to call the CylindricalVolumeBuilder
  /// call operator
  ///
  /// @note this function throws an exception if the provided internal
  /// blueprints do not fit into the restriction
  ///
  /// @return the new (cylindrical) container volume
  std::vector<std::shared_ptr<DetectorVolume>> operator()(
      const GeometryContext& gctx,
      const std::vector<InternalBlueprint>& iBlueprints,
      const Extent& restriction = Extent(false),
      const std::string& name = "Volume") noexcept(false);

  /// Private access method to the logger
  const Logger& logger() const { return *m_logger; }

  /// the logging instance
  std::unique_ptr<const Logger> m_logger;
};

struct CylindricalVolumeBuilderZ {
  /// Constructor with default logger
  ///
  /// @param level is the logging level
  CylindricalVolumeBuilderZ(Logging::Level level = Logging::INFO)
      : m_logger(getDefaultLogger("CylindricalVolumeBuilderZ", level)) {}

  /// Call operator that creates volumes for a
  /// cylindrical container ordered in Z
  ///
  /// @param gctx the geometry context for this building call
  /// @param iBlueprints the blue prints of the internal volume structure
  /// @param restriction is an extend that restricts the volume
  /// @param name the base name of the volumes
  ///
  /// In case no internal or only one internal blueprint is handed over
  /// this function will attempt to call the CylindricalVolumeBuilder
  /// call operator
  ///
  /// @note this function throws an exception if the provided internal
  /// blueprints do not fit into the restriction
  ///
  /// @return the new cylindrical volumes
  std::vector<std::shared_ptr<DetectorVolume>> operator()(
      const GeometryContext& gctx,
      const std::vector<InternalBlueprint>& iBlueprints,
      const Extent& restriction = Extent(false),
      const std::string& name = "Volume") noexcept(false);

  /// Private access method to the logger
  const Logger& logger() const { return *m_logger; }

  /// the logging instance
  std::unique_ptr<const Logger> m_logger;
};

struct CylindricalVolumeBuilderPhi {
  /// Constructor with default logger
  ///
  /// @param level is the logging level
  CylindricalVolumeBuilderPhi(Logging::Level level = Logging::INFO)
      : m_logger(getDefaultLogger("CylindricalVolumeBuilderPhi", level)) {}

  /// Call operator that creates volumes for a
  /// cylindrical container ordered in Phi
  ///
  /// @param gctx the geometry context for this building call
  /// @param iBlueprints the blue prints of the internal volume structure
  /// @param restriction is an extend that restricts the volume
  /// @param name the base name of the volumes
  ///
  /// In case no internal or only one internal blueprint is handed over
  /// this function will attempt to call the CylindricalVolumeBuilder
  /// call operator
  ///
  /// @note this function throws an exception if the provided internal
  /// blueprints do not fit into the restriction
  ///
  /// @return the new cylindrical volumes
  std::vector<std::shared_ptr<DetectorVolume>> operator()(
      const GeometryContext& gctx,
      const std::vector<InternalBlueprint>& iBlueprints,
      const Extent& restriction = Extent(false),
      const std::string& name = "Volume") noexcept(false);

  /// Private access method to the logger
  const Logger& logger() const { return *m_logger; }

  /// the logging instance
  std::unique_ptr<const Logger> m_logger;
};

}  // namespace Acts
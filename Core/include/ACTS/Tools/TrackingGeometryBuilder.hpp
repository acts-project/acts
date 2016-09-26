// This file is part of the ACTS project.
//
// Copyright (C) 2016 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// TrackingGeometryBuilder.h, ACTS project
///////////////////////////////////////////////////////////////////

#ifndef ACTS_TOOLS_CYLINDERGEOMETRYBUILDER_H
#define ACTS_TOOLS_CYLINDERGEOMETRYBUILDER_H 1

#include <list>
#include <memory>
#include "ACTS/Tools/ITrackingGeometryBuilder.hpp"
#include "ACTS/Tools/ITrackingVolumeBuilder.hpp"
#include "ACTS/Tools/ITrackingVolumeHelper.hpp"
#include "ACTS/Utilities/Definitions.hpp"
#include "ACTS/Utilities/Logger.hpp"

namespace Acts {
class TrackingGeometry;

/// @class GeometryBuilder
///
/// The Acts::TrackingGeometry Builder for volumes that wrap around another
///
/// It retrieves ITrackingVolumeBuilders as configured and builds one
/// detector around the other one.
///
class TrackingGeometryBuilder : public ITrackingGeometryBuilder
{
public:
  /// @struct Config
  /// Nested Configuration for the CylinderVolumeBuilder
  struct Config
  {
    /// the list of trackign volume builders
    std::list<std::shared_ptr<ITrackingVolumeBuilder>> trackingVolumeBuilders{};
    /// the tracking volume helper for detector construction
    std::shared_ptr<ITrackingVolumeHelper> trackingVolumeHelper = nullptr;
  };

  /// Constructor
  /// @param cgb is the configuration struct for this builder
  /// @param logger logging instance
  TrackingGeometryBuilder(const Config&           cgbConfig,
                          std::unique_ptr<Logger> logger
                          = getDefaultLogger("TrackingGeometryBuilder",
                                             Logging::INFO));

  /// Destructor
  virtual ~TrackingGeometryBuilder() = default;

  /// TrackingGeometry Interface method
  /// returns a unique pointer to a TrackingGeometry
  virtual std::unique_ptr<TrackingGeometry>
  trackingGeometry() const final;

  /// Set configuration method
  /// @param cgbConfig is the new configuration struct
  void
  setConfiguration(const Config& cgbConfig);

  /// Get configuration method
  /// @return the current configuration
  Config
  getConfiguration() const;

  /// set logging instance
  void
  setLogger(std::unique_ptr<Logger> logger);

private:
  /// Configuration member
  Config m_cfg;

  /// Private access method to the logger
  const Logger&
  logger() const
  {
    return *m_logger;
  }

  /// the logging instance
  std::unique_ptr<Logger> m_logger;
};

inline TrackingGeometryBuilder::Config
TrackingGeometryBuilder::getConfiguration() const
{
  return m_cfg;
}

}  // end of namespace

#endif  // ACTS_TOOLS_CYLINDERGEOMETRYBUILDER_H

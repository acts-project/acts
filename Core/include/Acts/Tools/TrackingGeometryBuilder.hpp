// This file is part of the Acts project.
//
// Copyright (C) 2016-2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// TrackingGeometryBuilder.h, Acts project
///////////////////////////////////////////////////////////////////

#pragma once
#include <list>
#include <memory>
#include "Acts/Tools/ITrackingGeometryBuilder.hpp"
#include "Acts/Tools/ITrackingVolumeBuilder.hpp"
#include "Acts/Tools/ITrackingVolumeHelper.hpp"
#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Utilities/Logger.hpp"

namespace Acts {
class TrackingGeometry;

/// @class GeometryBuilder
///
/// The Acts::TrackingGeometry Builder for volumes that wrap around another
///
/// It retrieves an array of ITrackingVolumeBuilder tools that are configured
/// to be built in sequence, where the output of one volume builder is provided
/// to the next volume volume builder and accordingly
/// - contained (e.g. a final insertion of a beam pipe of longer extend)
/// - wrapped (e.g. an outer detector wrapping an inner one)
/// - attached (e.g. a neighbor detector attaching to the previous one)
///
/// The returned volume of each step must be processable by the previous step
class TrackingGeometryBuilder : public ITrackingGeometryBuilder
{
public:
  /// @struct Config
  /// Nested Configuration for the CylinderVolumeBuilder
  struct Config
  {

    /// the list of tracking volume builders
    std::list<std::shared_ptr<const ITrackingVolumeBuilder>>
        trackingVolumeBuilders{};

    /// the tracking volume helper for detector construction
    std::shared_ptr<const ITrackingVolumeHelper> trackingVolumeHelper = nullptr;
  };

  /// Constructor
  ///
  /// @param [in] cgbConfig is the configuration struct for this builder
  /// @param [in] logger logging instance
  TrackingGeometryBuilder(const Config&                 cgbConfig,
                          std::unique_ptr<const Logger> logger
                          = getDefaultLogger("TrackingGeometryBuilder",
                                             Logging::INFO));

  /// Destructor
  ~TrackingGeometryBuilder() override = default;

  /// TrackingGeometry Interface method
  /// @return a unique pointer to a TrackingGeometry
  std::unique_ptr<const TrackingGeometry>
  trackingGeometry() const final;

  /// Set configuration method
  ///
  /// @param cgbConfig is the new configuration struct
  void
  setConfiguration(const Config& cgbConfig);

  /// Get configuration method
  /// @return the current configuration
  Config
  getConfiguration() const;

  /// set logging instance
  /// @param newLogger the new logging instance
  void
  setLogger(std::unique_ptr<const Logger> newLogger);

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
  std::unique_ptr<const Logger> m_logger;
};

inline TrackingGeometryBuilder::Config
TrackingGeometryBuilder::getConfiguration() const
{
  return m_cfg;
}

}  // namespace
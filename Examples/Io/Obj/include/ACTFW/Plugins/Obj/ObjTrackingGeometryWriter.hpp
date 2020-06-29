// This file is part of the Acts project.
//
// Copyright (C) 2017-2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <Acts/Utilities/Logger.hpp>
#include <Acts/Visualization/ViewConfig.hpp>
#include <fstream>
#include <iostream>
#include <mutex>

#include "ACTFW/Framework/AlgorithmContext.hpp"
#include "ACTFW/Framework/ProcessCode.hpp"

namespace Acts {
class TrackingVolume;
class TrackingGeometry;
}  // namespace Acts

namespace FW {
namespace Obj {

/// @class ObjTrackingGeometryWriter
///
/// An Obj writer for the geometry: TrackingGeometry master
/// It delegates the writing of surfaces to the surface writers
class ObjTrackingGeometryWriter {
 public:
  // @class Config
  //
  // The nested config class
  class Config {
   public:
    std::shared_ptr<const Acts::Logger> logger;

    std::string name = "";

    double outputScalor = 1.0;   ///< scale output values
    size_t outputPrecision = 6;  ///< floating point precision

    Acts::ViewConfig containerView = Acts::ViewConfig({220, 220, 220});
    Acts::ViewConfig volumeView = Acts::ViewConfig({220, 220, 0});
    Acts::ViewConfig sensitiveView = Acts::ViewConfig({0, 180, 240});
    Acts::ViewConfig passiveView = Acts::ViewConfig({240, 280, 0});
    Acts::ViewConfig gridView = Acts::ViewConfig({220, 0, 0});

    Config(const std::string& lname = "ObjTrackingGeometryWriter",
           Acts::Logging::Level lvl = Acts::Logging::INFO)
        : logger(Acts::getDefaultLogger(lname, lvl)), name(lname) {}
  };

  /// Constructor
  /// @param cfg is the configuration class
  ObjTrackingGeometryWriter(const Config& cfg);

  /// Framework name() method
  /// @return the name of the tool
  std::string name() const;

  /// The write interface
  /// @param context the Algorithm/Event context of this call
  /// @param tGeometry is the geometry to be written out
  /// @return ProcessCode to indicate success/failure
  FW::ProcessCode write(const AlgorithmContext& context,
                        const Acts::TrackingGeometry& tGeometry);

 private:
  Config m_cfg;  ///< the config class

  /// process this volume
  /// @param context the Algorithm/Event context for this call
  /// @param tVolume the volume to be processed
  void write(const AlgorithmContext& context,
             const Acts::TrackingVolume& tVolume);

  /// Private access to the logging instance
  const Acts::Logger& logger() const { return *m_cfg.logger; }
};

}  // namespace Obj
}  // namespace FW

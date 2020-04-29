// This file is part of the Acts project.
//
// Copyright (C) 2017-2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <Acts/Surfaces/Surface.hpp>
#include <Acts/Utilities/Logger.hpp>
#include <fstream>
#include <iostream>
#include <mutex>

#include "ACTFW/Framework/AlgorithmContext.hpp"
#include "ACTFW/Framework/ProcessCode.hpp"
#include "ACTFW/Plugins/Obj/ObjHelper.hpp"

namespace FW {
namespace Obj {

/// @class ObjSurfaceWriter
///
/// An Obj writer for the geometry: surface section
///
class ObjSurfaceWriter {
 public:
  // @class Config
  //
  // The nested config class for the Surface writer
  class Config {
   public:
    /// the default logger
    std::shared_ptr<const Acts::Logger> logger;
    /// the name of the algorithm
    std::string name;
    /// approximate cyinders by that
    unsigned int outputPhiSegemnts = 72;
    /// write thickness if available
    double outputThickness = 2.;
    /// write sensitive surfaces
    bool outputSensitive = true;
    /// write the layer surface out
    bool outputLayerSurface = true;
    /// output scalor
    double outputScalor = 1.;
    /// precision for out
    unsigned int outputPrecision = 6;
    /// file prefix to be written out
    std::string filePrefix = "";
    /// prefixes
    /// @todo These aren't used anywhere, should they be dropped?
    std::string planarPrefix = "";
    std::string cylinderPrefix = "";
    std::string diskPrefix = "";
    /// the output stream
    std::shared_ptr<std::ofstream> outputStream = nullptr;

    Config(const std::string& lname = "ObjSurfaceWriter",
           Acts::Logging::Level lvl = Acts::Logging::INFO)
        : logger(Acts::getDefaultLogger(lname, lvl)), name(lname) {}
  };

  /// Constructor
  ///
  /// @param cfg is the configuration class
  ObjSurfaceWriter(const Config& cfg);

  /// Framework name() method
  std::string name() const;

  /// The write interface
  /// @param context the Algorithm/Event context of this call
  /// @param surface to be written out
  FW::ProcessCode write(const AlgorithmContext& context,
                        const Acts::Surface& surface);

  /// write a bit of string
  /// @param is the string to be written
  FW::ProcessCode write(const std::string& sinfo);

 private:
  Config m_cfg;                  ///< the config class
  Obj::VtnCounter m_vtnCounter;  ///< vertex, texture, normal
  std::mutex m_write_mutex;      ///< mutex to protect multi-threaded writes

  /// Private access to the logging instance
  const Acts::Logger& logger() const { return *m_cfg.logger; }
};

inline FW::ProcessCode ObjSurfaceWriter::write(const std::string& sinfo) {
  // lock the mutex for writing
  std::lock_guard<std::mutex> lock(m_write_mutex);
  // and write
  (*m_cfg.outputStream) << sinfo;
  return FW::ProcessCode::SUCCESS;
}

}  // namespace Obj
}  // namespace FW

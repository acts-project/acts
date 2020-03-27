// This file is part of the Acts project.
//
// Copyright (C) 2017-2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// RootMaterialWriter.h
///////////////////////////////////////////////////////////////////

#pragma once

#include <Acts/Geometry/GeometryID.hpp>
#include <Acts/Geometry/TrackingGeometry.hpp>
#include <Acts/Geometry/TrackingVolume.hpp>
#include <Acts/Material/IMaterialDecorator.hpp>
#include <Acts/Material/ISurfaceMaterial.hpp>
#include <Acts/Material/IVolumeMaterial.hpp>
#include <Acts/Surfaces/Surface.hpp>
#include <Acts/Utilities/Definitions.hpp>
#include <Acts/Utilities/Logger.hpp>
#include <map>
#include <mutex>

#include "ACTFW/Framework/ProcessCode.hpp"

namespace Acts {
using SurfaceMaterialMap =
    std::map<GeometryID, std::shared_ptr<const ISurfaceMaterial>>;
using VolumeMaterialMap =
    std::map<GeometryID, std::shared_ptr<const IVolumeMaterial>>;
using DetectorMaterialMaps = std::pair<SurfaceMaterialMap, VolumeMaterialMap>;
}  // namespace Acts

namespace FW {

/// @brief Material decorator from Root format
///
/// This reads in material maps for surfaces and volumes
/// from a root file
class RootMaterialWriter {
 public:
  /// @class Config
  ///
  /// Configuration of the Writer
  struct Config {
    /// Steering to handle sensitive data
    bool processSensitives = true;

    /// Steering to handle approach data
    bool processApproaches = true;

    /// Steering to handle representing data
    bool processRepresenting = true;

    /// Steering to handle boundary data
    bool processBoundaries = true;

    /// Steering to handle volume data
    bool processVolumes = true;

    /// The name of the output tree
    std::string folderNameBase = "Material";
    /// The volume identification string
    std::string voltag = "_vol";
    /// The boundary identification string
    std::string boutag = "_bou";
    /// The layer identification string
    std::string laytag = "_lay";
    /// The approach identification string
    std::string apptag = "_app";
    /// The sensitive identification string
    std::string sentag = "_sen";
    /// The bin number tag
    std::string ntag = "n";
    /// The value tag -> binning values: binZ, binR, binPhi, etc.
    std::string vtag = "v";
    /// The option tag -> binning options: open, closed
    std::string otag = "o";
    /// The range min tag: min value
    std::string mintag = "min";
    /// The range max tag: max value
    std::string maxtag = "max";
    /// The thickness tag
    std::string ttag = "t";
    /// The x0 tag
    std::string x0tag = "x0";
    /// The l0 tag
    std::string l0tag = "l0";
    /// The A tag
    std::string atag = "A";
    /// The Z tag
    std::string ztag = "Z";
    /// The rho tag
    std::string rhotag = "rho";
    /// The name of the output file
    std::string fileName = "material-maps.root";
    /// The default logger
    std::shared_ptr<const Acts::Logger> logger;
    // The name of the writer
    std::string name = "";

    /// Constructor
    ///
    /// @param lname Name of the writer tool
    /// @param lvl The output logging level
    Config(const std::string& lname = "RootMaterialWriter",
           Acts::Logging::Level lvl = Acts::Logging::INFO)
        : logger(Acts::getDefaultLogger(lname, lvl)), name(lname) {}
  };

  /// Constructor
  ///
  /// @param cfg The configuration struct
  RootMaterialWriter(const Config& cfg);

  /// Virtual destructor
  ~RootMaterialWriter() = default;

  /// Write out the material map
  ///
  /// @param detMaterial is the SurfaceMaterial and VolumeMaterial maps
  void write(const Acts::DetectorMaterialMaps& detMaterial);

  /// Write out the material map from Geometry
  ///
  /// @param tGeometry is the TrackingGeometry
  void write(const Acts::TrackingGeometry& tGeometry);

 private:
  /// Collect the material from the tracking geometry
  ///
  /// @param tVolume The TrackingVolume for the material to be collected
  /// @param [in,out] detMatMap the map to be filled
  void collectMaterial(const Acts::TrackingVolume& tVolume,
                       Acts::DetectorMaterialMaps& detMatMap);

  /// Collect the material from the tracking geometry
  ///
  /// @param tLayer The TrackingVolume for the material to be collected
  /// @param [in,out] detMatMap the map to be filled
  void collectMaterial(const Acts::Layer& tLayer,
                       Acts::DetectorMaterialMaps& detMatMap);

  /// The config class
  Config m_cfg;

  /// Private access to the logging instance
  const Acts::Logger& logger() const { return *m_cfg.logger; }
};

}  // namespace FW

// This file is part of the Acts project.
//
// Copyright (C) 2017-2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <Acts/Geometry/GeometryID.hpp>
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

class TFile;

namespace Acts {
using SurfaceMaterialMap =
    std::map<GeometryID, std::shared_ptr<const ISurfaceMaterial>>;
using VolumeMaterialMap =
    std::map<GeometryID, std::shared_ptr<const IVolumeMaterial>>;
}  // namespace Acts

namespace FW {

/// @class RootMaterialDecorator
///
/// @brief Read the collection of SurfaceMaterial & VolumeMaterial
class RootMaterialDecorator : public Acts::IMaterialDecorator {
 public:
  /// @class Config
  /// Configuration of the Reader
  class Config {
   public:
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
    Config(const std::string& lname = "MaterialReader",
           Acts::Logging::Level lvl = Acts::Logging::INFO)
        : logger(Acts::getDefaultLogger(lname, lvl)), name(lname) {}
  };

  /// Constructor
  ///
  /// @param cfg configuration struct for the reader
  RootMaterialDecorator(const Config& cfg);

  /// Destructor
  ~RootMaterialDecorator();

  /// Decorate a surface
  ///
  /// @param surface the non-cost surface that is decorated
  void decorate(Acts::Surface& surface) const final {
    // Null out the material for this surface
    if (m_clearSurfaceMaterial) {
      surface.assignSurfaceMaterial(nullptr);
    }
    // Try to find the surface in the map
    auto sMaterial = m_surfaceMaterialMap.find(surface.geoID());
    if (sMaterial != m_surfaceMaterialMap.end()) {
      surface.assignSurfaceMaterial(sMaterial->second);
    }
  }

  /// Decorate a TrackingVolume
  ///
  /// @param volume the non-cost volume that is decorated
  void decorate(Acts::TrackingVolume& volume) const final {
    // Null out the material for this volume
    if (m_clearSurfaceMaterial) {
      volume.assignVolumeMaterial(nullptr);
    }
    // Try to find the surface in the map
    auto vMaterial = m_volumeMaterialMap.find(volume.geoID());
    if (vMaterial != m_volumeMaterialMap.end()) {
      volume.assignVolumeMaterial(vMaterial->second);
    }
  }

 private:
  /// The config class
  Config m_cfg;

  /// The input file
  TFile* m_inputFile{nullptr};

  /// Surface based material
  Acts::SurfaceMaterialMap m_surfaceMaterialMap;

  /// Volume based material
  Acts::VolumeMaterialMap m_volumeMaterialMap;

  bool m_clearSurfaceMaterial{true};

  /// Private access to the logging instance
  const Acts::Logger& logger() const { return *m_cfg.logger; }
};

}  // namespace FW

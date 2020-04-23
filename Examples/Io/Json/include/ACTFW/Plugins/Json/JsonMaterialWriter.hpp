// This file is part of the Acts project.
//
// Copyright (C) 2017-2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// JsonMaterialWriter.h
///////////////////////////////////////////////////////////////////

#pragma once

#include <mutex>

#include "ACTFW/Framework/ProcessCode.hpp"
#include "Acts/Geometry/GeometryID.hpp"
#include "Acts/Material/ISurfaceMaterial.hpp"
#include "Acts/Material/IVolumeMaterial.hpp"
#include "Acts/Plugins/Json/JsonGeometryConverter.hpp"
#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Utilities/Logger.hpp"

namespace Acts {

class TrackingGeometry;

using SurfaceMaterialMap =
    std::map<GeometryID, std::shared_ptr<const ISurfaceMaterial>>;

using VolumeMaterialMap =
    std::map<GeometryID, std::shared_ptr<const IVolumeMaterial>>;

using DetectorMaterialMaps = std::pair<SurfaceMaterialMap, VolumeMaterialMap>;
}  // namespace Acts

namespace FW {

namespace Json {

/// @class Json Material writer
///
/// @brief Writes out Detector material maps
/// using the Json Geometry converter
class JsonMaterialWriter {
 public:
  /// Constructor
  ///
  /// @param cfg The configuration struct of the converter
  JsonMaterialWriter(const Acts::JsonGeometryConverter::Config& cfg,
                     const std::string& fileName);

  /// Virtual destructor
  ~JsonMaterialWriter();

  /// Write out the material map
  ///
  /// @param detMaterial is the SurfaceMaterial and VolumeMaterial maps
  void write(const Acts::DetectorMaterialMaps& detMaterial);

  /// Write out the material map from Geometry
  ///
  /// @param tGeometry is the TrackingGeometry
  void write(const Acts::TrackingGeometry& tGeometry);

 private:
  /// The config class of the converter
  Acts::JsonGeometryConverter::Config m_cfg;

  /// The file name
  std::string m_fileName;

  /// Private access to the logging instance
  const Acts::Logger& logger() const { return *m_cfg.logger; }
};

}  // namespace Json
}  // namespace FW

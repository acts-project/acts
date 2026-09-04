// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Material/IMaterialDecorator.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/MaterialMapping/IMaterialWriter.hpp"
#include "ActsPlugins/Root/RootMaterialMapIo.hpp"

#include <memory>
#include <string>

class TFile;

namespace Acts {
class ISurfaceMaterial;
class IVolumeMaterial;
class Layer;
class TrackingGeometry;
class TrackingVolume;

}  // namespace Acts

namespace ActsExamples {

/// @brief Material decorator from Root format
///
/// This reads in material maps for surfaces and volumes
/// from a root file
class RootMaterialWriter : public IMaterialWriter {
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

    /// The accessor configuration
    ActsPlugins::RootMaterialMapIo::Config accessorConfig;

    /// The accessor options
    ActsPlugins::RootMaterialMapIo::Options accessorOptions;

    /// The name of the output file
    std::string filePath = "material-maps.root";
    /// The file mode
    std::string fileMode = "RECREATE";
  };

  /// Constructor
  ///
  /// @param config The configuration struct
  /// @param level The log level
  RootMaterialWriter(const Config& config, Acts::Logging::Level level);

  /// Virtual destructor
  ~RootMaterialWriter() override;

  /// Write out the material map
  ///
  /// @param detMaterial is the SurfaceMaterial and VolumeMaterial maps
  void writeMaterial(
      const Acts::TrackingGeometryMaterial& detMaterial) override;

  /// Write out the material map from Geometry
  ///
  /// @param tGeometry is the TrackingGeometry
  void write(const Acts::TrackingGeometry& tGeometry);

  /// Get readonly access to the config parameters
  const Config& config() const { return m_cfg; }

 private:
  /// Collect the material from the tracking geometry
  ///
  /// @param tVolume The TrackingVolume for the material to be collected
  /// @param [in,out] detMatMap the map to be filled
  void collectMaterial(const Acts::TrackingVolume& tVolume,
                       Acts::TrackingGeometryMaterial& detMatMap);

  /// Collect the material from the tracking geometry
  ///
  /// @param tLayer The TrackingVolume for the material to be collected
  /// @param [in,out] detMatMap the map to be filled
  void collectMaterial(const Acts::Layer& tLayer,
                       Acts::TrackingGeometryMaterial& detMatMap);

  /// The config class
  Config m_cfg;

  /// The logger instance
  std::unique_ptr<const Acts::Logger> m_logger;

  /// Private access to the logging instance
  const Acts::Logger& logger() const { return *m_logger; }

  TFile* m_outputFile{nullptr};
};

}  // namespace ActsExamples

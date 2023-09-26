// This file is part of the Acts project.
//
// Copyright (C) 2022-2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Detector/DetectorVolume.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Detector/DetectorComponents.hpp"
#include "Acts/Detector/ProtoBinning.hpp"
#include "Acts/Detector/DetectorComponents.hpp"
#include "Acts/Utilities/Logger.hpp"

#include <iostream>
#include <memory>
#include <string>
#include <vector>
#include <optional>

class G4VPhysicalVolume;

namespace Acts {
namespace Experimental {
class DetectorVolume;
}  // namespace Experimental
}  // namespace Acts

namespace ActsExamples {

///@class MockupsectorBuilder
///
/// This is a class that builds a chamber by reading a gdml file

class MockupSectorBuilder {
 public:
  /// Nested configuration struct
  struct Config {
    // The path of the gdml file that holds the mockup geometry
    std::string gdmlPath = "";

    // The number of sectors we want to build 
    int NumberOfSectors = 1;

    float toleranceOverlap = 10.;

   // Flag for the multilayer binning 
    bool robustMode = false;

  };

  /// Nested configuration struct for a multi layer of a muon spectrometer chamber
  struct MultiLayerConfig {
    
    // The name of the multi layer
    std::string name;

    // The names of the sensitive surfaces
    std::vector<std::string> SensitiveNames;

    // The names of the passive surfaces
    std::vector<std::string> PassiveNames;   

  };

  /// Nested configuration for chamber construction
  struct ChamberConfig{

    //The name of the chamber
    std::string name;

    //The internal volumes of the chamber (e.g the two multi-layer detector volumes for a muon spectrometer chamber)
    std::vector<std::shared_ptr<Acts::Experimental::DetectorVolume>> internalVolumes;

  };


  /// Constructor
  ///@param config The configuration struct
  MockupSectorBuilder(const Config& config, std::unique_ptr<const Acts::Logger> logger =
          Acts::getDefaultLogger("MockupSectorBuilder", Acts::Logging::VERBOSE));

  /// Destructor
  ~MockupSectorBuilder() = default;

  /// Build MultiLayer
  /// @param mlConfig The configuration for a multilayer volume
  std::shared_ptr<Acts::Experimental::DetectorComponent> buildMultiLayer(
       MultiLayerConfig& mlConfig);

  /// Build chamber
  /// @param detComp The Detector Components to build the chamber (e.g the multi layers)
  std::shared_ptr<Acts::Experimental::DetectorComponent> buildChamber(
    const ChamberConfig& chConfig);

  /// Build Sector
  /// @param det_volumes The vector that contains the detector volumes of the Sector (e.g the chambers)
  std::shared_ptr<const Acts::Experimental::DetectorComponent> buildSector(
      std::vector<std::shared_ptr<Acts::Experimental::DetectorVolume>>
          detVolumes);

  /// Draw the sector in an obj file
  /// @param nameObjFile The name of the obj file where the sector will be saved
  void drawSector(const std::shared_ptr<Acts::Experimental::DetectorVolume>&
                      detectorVolumeSector,
                  const std::string& nameObjFile);

 private:

  ///Configuration member
  Config mCfg;

  /// Geant4 World
  G4VPhysicalVolume* g4World = nullptr;

  /// Maximum number of sctors to build
  int maxNumberOfSectors = 8;

  /// A binning definition function for binning along y and z axis
  /// @param binEdgesY The edges (min and max) for y axis
  /// @param binEdgesZ The edges(min and max) for z axis
  /// @param radius The radius of the tube surface for the bin width definition
  std::vector<Acts::Experimental::ProtoBinning> defineBinning(std::pair<double,double> binEdgesZ, 
    std::pair<double,double> binEdgesY, float radius);

  /// logging instance
  std::unique_ptr<const Acts::Logger> mLogger;

   /// Private access method to the logger
  const Acts::Logger& logger() const { return *mLogger; }

};

}  // namespace ActsExamples

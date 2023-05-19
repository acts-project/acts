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
#include "Acts/Navigation/SurfaceCandidatesUpdators.hpp"

#include <iostream>
#include <memory>
#include <string>
#include <vector>

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
  friend struct ChamberConfig;

 public:
  // Nested configuration struct
  struct Config {
    friend struct ChamberConfig;

    // The path of the gdml file that holds the mockup geometry
    std::string gdmlPath = "";

    // The number of sectors we want to create
    int NumberOfSectors = 1;

    float toleranceOverlap = 10.;
  };

  /// Nested configuration struct for chamber
  struct ChamberConfig {
    friend class MockupSectorBuilder;
    friend struct Config;

    // The name of the chamber
    std::string name;

    // The names of the sensitive surfaces
    std::vector<std::string> SensitiveNames;

    // The names of the passive surfaces
    std::vector<std::string> PassiveNames;

    // The indexed surfaces
    // SurfaceCandidatesUpdator
    // getIndexedStraws(std::vector<std::shared_ptr<Acts::Surface>>);

    // Function that returns the straws of the chamber of a specific chamber
    std::vector<std::shared_ptr<Acts::Surface>> getStraws() {
      return ChambersStraws;
    };

   private:
    // void setStraws
    std::vector<std::shared_ptr<Acts::Surface>> ChambersStraws;
  };

  /// Constructor
  ///@param config The configuration struct
  MockupSectorBuilder(const Config& config);

  /// Destructor
  ~MockupSectorBuilder() = default;

  /// Build chamber
  /// @param gctx The current geometry context object
  /// @param chamber_config The configuration chamber struct
  std::shared_ptr<Acts::Experimental::DetectorVolume> buildChamber(
      ChamberConfig& chamberConfig);

  /// Build Sector
  /// @param det_volumes The vector that contains the detector volumes of the Sector
  /// @param gctx The current geometry context object
  std::shared_ptr<Acts::Experimental::DetectorVolume> buildSector(
      std::vector<std::shared_ptr<Acts::Experimental::DetectorVolume>>
          detVolumes);

  /// Draw the sector in an obj file
  /// @param nameObjFile The name of the obj file where the sector will be saved
  void drawSector(const std::shared_ptr<Acts::Experimental::DetectorVolume>&
                      detectorVolumeSector,
                  const std::string& nameObjFile);

 private:
  Config mCfg;

  G4VPhysicalVolume* g4World = nullptr;

  int maxNumberOfSectors = 8;

  /// Build multi-layer
  /// @param strawSurfaces The straw surfaces from the chamber
  /// @param chamber_pos The position of the chamber that includes the multilayers
  /// @param multilayer_name The name of the multilayer (wether it is the up or the bottom multilayer)
  std::shared_ptr<Acts::Experimental::DetectorVolume> buildMultiLayer(
      const std::vector<std::shared_ptr<Acts::Surface>> strawSurfaces,
      std::string multilayer_name);
};

}  // namespace ActsExamples

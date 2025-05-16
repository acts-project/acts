// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Plugins/GeoModel/GeoModelTree.hpp"
#include "ActsExamples/Geant4/Geant4ConstructionOptions.hpp"
#include "ActsExamples/Geant4/SensitiveSurfaceMapper.hpp"
#include <G4VUserDetectorConstruction.hh>

#include <unordered_set>

class G4VPhysicalVolume;

namespace ActsExamples {

/// Construct the Geant4 detector from a GeoModel world volume
class GeoModelGeant4DetectorConstruction final
    : public G4VUserDetectorConstruction {
 public:
  /** @brief configuration object to steer the GeoModel -> G4 geometry translation. */
  struct Config: public Geant4ConstructionOptions{
      Config() = default;

      Config(const Geant4ConstructionOptions& opts):
          Geant4ConstructionOptions{opts}{}
      
      std::vector<std::string> sensitiveVols{};

      Acts::GeoModelTree geoModelTree{};

  };
 
  /// @param geoModelTree is the GeoModel tree containing the world volume
  /// @param regionCreators are the region creators
  GeoModelGeant4DetectorConstruction(const Config& cfg);

  /// Read the file and parse it to construct the Geant4 description
  ///
  /// @note to simplify further setup withiin the ACTS framework
  /// the world is cached when first constructed
  G4VPhysicalVolume* Construct() override;

 private:
  /** @brief Prefix to be added to the volume name string in order to later record the G4 volume
   *         by the ActsExamples::geant4::SensitiveSteppingAction */
  constexpr static std::string_view mappingPrefix = ActsExamples::Geant4::SensitiveSurfaceMapper::mappingPrefix;
  /** @brief Check whether the volume's name satisfies at least one of the sensitive volume
   *         tokens and then prepends the mapping prefix. Otherwise, the test is executed 
   *         on all daughter volumes. As the GeoModel -> G4 translation also respects the 
   *         volume sharing, the processed set keeps track whether this particular instance
   *         has already been tested.
   *  @param g4Vol: Pointer to the G4 logical volume to be tested for being sensitive
   *  @param processed: Reference to the container keeping track of the tested volumes */
  void markSensitiveVols(G4VPhysicalVolume* g4Vol,
                         std::unordered_set<G4LogicalVolume*>& processed) const;
  /// Construction options
  Config m_options;
  /// The world volume
  G4VPhysicalVolume* m_g4World = nullptr;
};

}  // namespace ActsExamples

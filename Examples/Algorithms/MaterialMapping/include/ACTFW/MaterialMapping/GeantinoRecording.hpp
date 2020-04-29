// This file is part of the Acts project.
//
// Copyright (C) 2017 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <Acts/Propagator/MaterialInteractor.hpp>
#include <Acts/Utilities/Definitions.hpp>
#include <Acts/Utilities/Logger.hpp>
#include <G4RunManager.hh>
#include <memory>

#include "ACTFW/Framework/BareAlgorithm.hpp"
#include "ACTFW/Plugins/DD4hepG4/DD4hepToG4Svc.hpp"

namespace Acts {
// Using some short hands for Recorded Material
using RecordedMaterial = MaterialInteractor::result_type;

// And recorded material track
// - this is start:  position, start momentum
//   and the Recorded material
using RecordedMaterialTrack =
    std::pair<std::pair<Acts::Vector3D, Acts::Vector3D>, RecordedMaterial>;
}  // namespace Acts

namespace FW {

/// @class GeantinoRecording
///
/// @brief records the simulation geometry by using geantinos
///
/// The GeantinoRecording algorithm initiates the geant4 simulation
/// of the Geant4materialMapping plugin, creates and writes out
/// the MaterialTrack entities which are needed for material mapping.
/// The input for the geant4 geometry can be either provided by a gdml file
/// or an implementation of the IGeant4Service.
class GeantinoRecording : public FW::BareAlgorithm {
 public:
  /// @class Config
  struct Config {
    std::string geantMaterialCollection = "geant-material-tracks";

    /// The service possibly providing the Geant4 geometry (optional)
    /// @note If this is not set, the geometry should be given by gdml file
    std::shared_ptr<DD4hepG4::DD4hepToG4Svc> geant4Service = nullptr;

    /// The possible gmdl input (optional)
    std::string gdmlFile;
    /// The number of tracks per event
    size_t tracksPerEvent = 0;

    /// random number seed 1
    int seed1 = 12345;
    /// random number seed 2
    int seed2 = 45678;
  };

  /// Constructor
  GeantinoRecording(const Config& cnf,
                    Acts::Logging::Level level = Acts::Logging::INFO);

  FW::ProcessCode execute(const AlgorithmContext& context) const final override;

 private:
  /// The config object
  Config m_cfg;
  /// G4 run manager
  std::unique_ptr<G4RunManager> m_runManager;
};
}  // namespace FW

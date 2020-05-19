// This file is part of the Acts project.
//
// Copyright (C) 2017-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <memory>
#include <mutex>

#include "ACTFW/Framework/BareAlgorithm.hpp"
#include "Acts/Propagator/MaterialInteractor.hpp"
#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Utilities/Logger.hpp"

class G4RunManager;
class G4VUserDetectorConstruction;

namespace ActsExamples {

using RecordedMaterial = Acts::MaterialInteractor::result_type;
/// A material track with start position and momentum and recorded material.
using RecordedMaterialTrack =
    std::pair<std::pair<Acts::Vector3D, Acts::Vector3D>, RecordedMaterial>;

/// Records the simulation geometry using geantinos.
///
/// This initiates the Geant4 simulation, and creates and writes out
/// the MaterialTrack entities which are needed for material mapping.
class GeantinoRecording final : public FW::BareAlgorithm {
 public:
  struct Config {
    /// Output collection for the generated material tracks.
    std::string outputMaterialTracks = "geant-material-tracks";
    /// Detector construction object.
    std::unique_ptr<G4VUserDetectorConstruction> detectorConstruction;
    /// The number of tracks per event.
    size_t tracksPerEvent = 0;
    /// random number seed 1.
    int seed1 = 12345;
    /// random number seed 2.
    int seed2 = 45678;
  };

  GeantinoRecording(Config&& cfg, Acts::Logging::Level lvl);
  ~GeantinoRecording();

  FW::ProcessCode execute(const FW::AlgorithmContext& ctx) const final override;

 private:
  Config m_cfg;
  std::unique_ptr<G4RunManager> m_runManager;
  // has to be mutable; algorithm interface enforces object constness
  mutable std::mutex m_runManagerLock;
};

}  // namespace ActsExamples

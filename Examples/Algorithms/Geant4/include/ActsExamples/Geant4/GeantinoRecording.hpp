// This file is part of the Acts project.
//
// Copyright (C) 2017-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Propagator/MaterialInteractor.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "Acts/Utilities/PolymorphicValue.hpp"
#include "ActsExamples/Framework/BareAlgorithm.hpp"
#include "ActsExamples/Framework/ProcessCode.hpp"
#include "ActsExamples/Geant4/G4DetectorConstructionFactory.hpp"
#include "ActsExamples/Geant4/PrimaryGeneratorAction.hpp"

#include <memory>
#include <mutex>
#include <string>

#include "G4VUserDetectorConstruction.hh"

class G4RunManager;

namespace ActsExamples {

using RecordedMaterial = Acts::MaterialInteractor::result_type;
/// A material track with start position and momentum and recorded material.
using RecordedMaterialTrack =
    std::pair<std::pair<Acts::Vector3, Acts::Vector3>, RecordedMaterial>;

/// Records the simulation geometry using geantinos.
///
/// This initiates the Geant4 simulation, and creates and writes out
/// the MaterialTrack entities which are needed for material mapping.
class GeantinoRecording final : public BareAlgorithm {
 public:
  struct Config {
    /// Output collection for the generated material tracks.
    std::string outputMaterialTracks = "geant-material-tracks";
    /// Detector construction object.
    std::shared_ptr<G4DetectorConstructionFactory> detectorConstructionFactory;
    /// The number of tracks per event.
    size_t tracksPerEvent = 0;
    /// Configuration of the generator action
    Geant4::PrimaryGeneratorAction::Config generationConfig;
  };

  GeantinoRecording(Config config, Acts::Logging::Level level);
  ~GeantinoRecording();

  ActsExamples::ProcessCode execute(
      const ActsExamples::AlgorithmContext& ctx) const final override;

  /// Readonly access to the configuration
  const Config& config() const { return m_cfg; }

 private:
  Config m_cfg;
  std::unique_ptr<G4RunManager> m_runManager;
  // has to be mutable; algorithm interface enforces object constness
  mutable std::mutex m_runManagerLock;
};

}  // namespace ActsExamples

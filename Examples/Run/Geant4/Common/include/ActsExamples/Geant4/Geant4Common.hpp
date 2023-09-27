// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "ActsExamples/Framework/Sequencer.hpp"
#include "ActsExamples/Utilities/OptionsFwd.hpp"

#include <memory>
#include <vector>

namespace Acts {
class TrackingGeometry;
class MagneticFieldProvider;
}  // namespace Acts

class G4VUserDetectorConstruction;
class G4RunManager;
class G4UserRunAction;
class G4UserEventAction;
class G4UserTrackingAction;
class G4UserSteppingAction;
class G4VUserPhysicsList;

namespace ActsExamples {

/// Main function for running Geant4 with a specific detector.
///
/// @param vars the parsed variables
/// @param sequencer the event sequencer
/// @param runManager the Geant4 runmanager to use
/// @param detector the detector to be used
/// @param runActions the list of Geant4 user run action
/// @param eventActions the list of Geant4 user event action
/// @param trackingActions the list of Geant4 user tracking action
/// @param steppingActions the list of Geant4 user stepping action
/// @param trackingGeometry the tracking geometry for the sennsitive mapping
/// @param magneticField the ACTS magnetic field to be wrapped
/// @param materialRecording boolean flag to run material mapping
///
void setupGeant4Simulation(
    const ActsExamples::Options::Variables& vars,
    ActsExamples::Sequencer& sequencer,
    std::shared_ptr<G4RunManager> runManager,
    std::unique_ptr<G4VUserDetectorConstruction> detector,
    std::vector<G4UserRunAction*> runActions = {},
    std::vector<G4UserEventAction*> eventActions = {},
    std::vector<G4UserTrackingAction*> trackingActions = {},
    std::vector<G4UserSteppingAction*> steppingActions = {},
    const std::shared_ptr<const Acts::TrackingGeometry>& trackingGeometry =
        nullptr,
    const std::shared_ptr<const Acts::MagneticFieldProvider>& magneticField =
        nullptr,
    bool materialRecording = false);

/// Specific setup: Material Recording
///
/// @param vars the parsed variables
/// @param sequencer the event sequencer
/// @param g4DetectorFactory is the detector to be used
int runMaterialRecording(const ActsExamples::Options::Variables& vars,
                         std::unique_ptr<G4VUserDetectorConstruction> detector);

/// Specific setup: Geant4 Simulation
///
/// @param vars the parsed variables
/// @param sequencer the event sequencer
/// @param g4DetectorFactory is the detector to be used
/// @param trackingGeometry the tracking geometry for the sennsitive mapping
int runGeant4Simulation(
    const ActsExamples::Options::Variables& vars,
    std::unique_ptr<G4VUserDetectorConstruction> detector,
    const std::shared_ptr<const Acts::TrackingGeometry>& trackingGeometry);

}  // namespace ActsExamples

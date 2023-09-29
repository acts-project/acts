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

class DetectorConstructionFactory;

/// Main function for running Geant4 material recording with a specific
/// detector.
///
/// @param vars the parsed variables
/// @param sequencer the event sequencer
/// @param detectorConstructionFactory is the detector construction factory to be used
/// @param trackingGeometry the tracking geometry for the sennsitive mapping
/// @param magneticField the ACTS magnetic field to be wrapped
///
void setupMaterialRecording(
    const ActsExamples::Options::Variables& vars,
    ActsExamples::Sequencer& sequencer,
    std::shared_ptr<DetectorConstructionFactory> detectorConstructionFactory);

/// Main function for running Geant4 simulation with a specific detector.
///
/// @param vars the parsed variables
/// @param sequencer the event sequencer
/// @param detectorConstructionFactory is the detector construction factory to be used
/// @param trackingGeometry the tracking geometry for the sennsitive mapping
/// @param magneticField the ACTS magnetic field to be wrapped
///
void setupGeant4Simulation(
    const ActsExamples::Options::Variables& vars,
    ActsExamples::Sequencer& sequencer,
    std::shared_ptr<DetectorConstructionFactory> detectorConstructionFactory,
    std::shared_ptr<const Acts::TrackingGeometry> trackingGeometry,
    std::shared_ptr<const Acts::MagneticFieldProvider> magneticField);

/// Specific setup: Material Recording
///
/// @param vars the parsed variables
/// @param detectorConstructionFactory is the detector construction factory to be used
int runMaterialRecording(
    const ActsExamples::Options::Variables& vars,
    std::shared_ptr<ActsExamples::DetectorConstructionFactory>
        detectorConstructionFactory);

/// Specific setup: Geant4 Simulation
///
/// @param vars the parsed variables
/// @param detectorConstructionFactory is the detector construction factory to be used
/// @param trackingGeometry the tracking geometry for the sennsitive mapping
int runGeant4Simulation(
    const ActsExamples::Options::Variables& vars,
    std::shared_ptr<ActsExamples::DetectorConstructionFactory>
        detectorConstructionFactory,
    std::shared_ptr<const Acts::TrackingGeometry> trackingGeometry);

}  // namespace ActsExamples

// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Geometry/Detector.hpp"
#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/MagneticField/MagneticFieldProvider.hpp"
#include "Acts/Plugins/Geant4/Geant4DetectorElement.hpp"
#include "Acts/Plugins/Geant4/Geant4DetectorSurfaceFactory.hpp"
#include "Acts/Plugins/Geant4/Geant4PhysicalVolumeSelectors.hpp"
#include "Acts/Plugins/Python/Utilities.hpp"
#include "ActsExamples/Framework/IContextDecorator.hpp"
#include "ActsExamples/Geant4/GdmlDetectorConstruction.hpp"
#include "ActsExamples/Geant4/Geant4Simulation.hpp"
#include "ActsExamples/Geant4/MagneticFieldWrapper.hpp"
#include "ActsExamples/Geant4/MaterialPhysicsList.hpp"
#include "ActsExamples/Geant4/MaterialSteppingAction.hpp"
#include "ActsExamples/Geant4/ParticleTrackingAction.hpp"
#include "ActsExamples/Geant4/SensitiveSteppingAction.hpp"
#include "ActsExamples/Geant4/SensitiveSurfaceMapper.hpp"
#include "ActsExamples/Geant4/SimParticleTranslation.hpp"
#include "ActsExamples/Geant4Detector/Geant4Detector.hpp"
#include "ActsExamples/TelescopeDetector/TelescopeG4DetectorConstruction.hpp"

#include <memory>

#include <FTFP_BERT.hh>
#include <G4MagneticField.hh>
#include <G4RunManager.hh>
#include <G4UserEventAction.hh>
#include <G4UserRunAction.hh>
#include <G4UserSteppingAction.hh>
#include <G4UserTrackingAction.hh>
#include <G4VPhysicalVolume.hh>
#include <G4VUserDetectorConstruction.hh>
#include <G4VUserPrimaryGeneratorAction.hh>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;
using namespace pybind11::literals;

using namespace ActsExamples;
using namespace Acts;

namespace Acts::Python {
void addGeant4HepMC3(Context& ctx);
}

PYBIND11_MODULE(ActsPythonBindingsGeant4, mod) {
  py::class_<G4VUserDetectorConstruction>(mod, "G4VUserDetectorConstruction");

  py::class_<G4VPhysicalVolume>(mod, "G4VPhysicalVolume");

  // This is the actual class we're binding
  py::class_<GdmlDetectorConstruction, G4VUserDetectorConstruction>(
      mod, "GdmlDetectorConstructionImpl")
      .def("Construct", &GdmlDetectorConstruction::Construct);

  // This is a python-only factory method that returns the above class.
  // We can apply a return value policy here so that python does NOT assume
  // ownership of the returned pointer, and it is safe to pass to G4
  mod.def(
      "GdmlDetectorConstruction",
      [](const std::string& path) {
        return new GdmlDetectorConstruction(path);
      },
      py::return_value_policy::reference);

  py::class_<SensitiveSurfaceMapper, std::shared_ptr<SensitiveSurfaceMapper>>(
      mod, "SensitiveSurfaceMapper");

  ACTS_PYTHON_DECLARE_ALGORITHM(
      Geant4Simulation, mod, "Geant4Simulation", outputSimHits,
      outputParticlesInitial, outputParticlesFinal, outputMaterialTracks,
      randomNumbers, runManager, primaryGeneratorAction, runActions,
      eventActions, trackingActions, steppingActions, detectorConstruction,
      magneticField, sensitiveSurfaceMapper);

  mod.def(
      "materialRecordingConfig",
      [](Acts::Logging::Level level, G4VUserDetectorConstruction* detector,
         std::shared_ptr<const ActsExamples::RandomNumbers> randomNumbers,
         const std::string& inputParticles,
         const std::string& outputMaterialTracks) {
        // The Geant4 actions needed
        std::vector<G4UserRunAction*> runActions = {};
        std::vector<G4UserEventAction*> eventActions = {};
        std::vector<G4UserTrackingAction*> trackingActions = {};

        // Set the main Geant4 algorithm, primary generation, detector
        // construction
        Geant4Simulation::Config g4Cfg;
        g4Cfg.randomNumbers = std::move(randomNumbers);
        g4Cfg.runManager = std::make_shared<G4RunManager>();
        g4Cfg.runManager->SetUserInitialization(new MaterialPhysicsList(
            Acts::getDefaultLogger("MaterialPhysicsList", level)));

        MaterialSteppingAction::Config mStepCfg;
        mStepCfg.excludeMaterials = {"Air", "Vacuum"};
        std::vector<G4UserSteppingAction*> steppingActions = {
            new MaterialSteppingAction(
                mStepCfg,
                Acts::getDefaultLogger("MaterialSteppingAction", level))};

        // Read the particle from the generator
        SimParticleTranslation::Config g4PrCfg;
        g4PrCfg.inputParticles = inputParticles;
        g4PrCfg.forceParticle = true;
        g4PrCfg.forcedMass = 0.;
        g4PrCfg.forcedPdgCode = 999;
        // Set the material tracks at output
        g4Cfg.outputMaterialTracks = outputMaterialTracks;

        // Set the primarty generator
        g4Cfg.primaryGeneratorAction = new SimParticleTranslation(
            g4PrCfg, Acts::getDefaultLogger("SimParticleTranslation", level));
        g4Cfg.detectorConstruction = detector;

        // Set the user actions
        g4Cfg.runActions = runActions;
        g4Cfg.eventActions = eventActions;
        g4Cfg.trackingActions = trackingActions;
        g4Cfg.steppingActions = steppingActions;

        return g4Cfg;
      },
      "level"_a, "detector"_a, "randomNumbers"_a, "inputParticles"_a,
      "outputMaterialTracks"_a);

  mod.def(
      "geant4SimulationConfig",
      [](Acts::Logging::Level& level, G4VUserDetectorConstruction* detector,
         const std::string& inputParticles,
         const std::shared_ptr<const Acts::TrackingGeometry>& trackingGeometry,
         const std::shared_ptr<const Acts::MagneticFieldProvider>&
             magneticField,
         const std::vector<std::string>& volumeMappings,
         const std::vector<std::string>& materialMappings) {
        // The Geant4 actions needed
        std::vector<G4UserRunAction*> runActions = {};
        std::vector<G4UserEventAction*> eventActions = {};
        std::vector<G4UserTrackingAction*> trackingActions = {};
        std::vector<G4UserSteppingAction*> steppingActions = {};

        // Set the main Geant4 algorithm, primary generation, detector
        // construction
        Geant4Simulation::Config g4Cfg;

        g4Cfg.runManager = std::make_shared<G4RunManager>();
        g4Cfg.runManager->SetUserInitialization(new FTFP_BERT());

        ParticleTrackingAction::Config g4TrackCfg;
        ParticleTrackingAction* particleAction = new ParticleTrackingAction(
            g4TrackCfg,
            Acts::getDefaultLogger("ParticleTrackingAction", level));
        trackingActions.push_back(particleAction);

        SensitiveSteppingAction::Config g4StepCfg;
        SensitiveSteppingAction* sensitiveStepping =
            new SensitiveSteppingAction(
                g4StepCfg,
                Acts::getDefaultLogger("SensitiveSteppingAction", level));
        steppingActions.push_back(sensitiveStepping);

        // Read the particle from the generator
        SimParticleTranslation::Config g4PrCfg;
        g4PrCfg.inputParticles = inputParticles;

        // Set the primarty generator
        g4Cfg.primaryGeneratorAction = new SimParticleTranslation(
            g4PrCfg, Acts::getDefaultLogger("SimParticleTranslation", level));
        g4Cfg.detectorConstruction = detector;

        // Set the user actions
        g4Cfg.runActions = runActions;
        g4Cfg.eventActions = eventActions;
        g4Cfg.trackingActions = trackingActions;
        g4Cfg.steppingActions = steppingActions;

        // An ACTS Magnetic field is provided
        if (magneticField) {
          MagneticFieldWrapper::Config g4FieldCfg;
          g4FieldCfg.magneticField = magneticField;
          g4Cfg.magneticField = new MagneticFieldWrapper(g4FieldCfg);
        }

        // An ACTS TrackingGeometry is provided, so simulation for sensitive
        // detectors is turned on - they need to get matched first
        if (trackingGeometry) {
          SensitiveSurfaceMapper::Config ssmCfg;
          ssmCfg.trackingGeometry = trackingGeometry;

          // Take the default args if nothing provided
          if (not volumeMappings.empty()) {
            ssmCfg.volumeMappings = volumeMappings;
          }
          if (not materialMappings.empty()) {
            ssmCfg.materialMappings = materialMappings;
          }

          g4Cfg.sensitiveSurfaceMapper =
              std::make_shared<const SensitiveSurfaceMapper>(
                  ssmCfg,
                  Acts::getDefaultLogger("SensitiveSurfaceMapper", level));
        }

        return g4Cfg;
      },
      "level"_a, "detector"_a, "inputParticles"_a,
      py::arg("trackingGeometry") = nullptr, py::arg("magneticField") = nullptr,
      py::arg("volumeMappings") = std::vector<std::string>{},
      py::arg("materialMappings") = std::vector<std::string>{});

  {
    using Detector = ActsExamples::Telescope::TelescopeDetector;
    using DetectorConstruction =
        ActsExamples::Telescope::TelescopeG4DetectorConstruction;

    py::class_<DetectorConstruction, G4VUserDetectorConstruction>(
        mod, "TelescopeG4DetectorConstructionImpl");

    mod.def(
        "TelescopeG4DetectorConstruction",
        [](Detector& detector) {
          return new DetectorConstruction(detector.config);
        },
        py::return_value_policy::reference);
  }
  {
    using ISelector = Acts::IGeant4PhysicalVolumeSelector;
    auto is = py::class_<ISelector, std::shared_ptr<ISelector>>(
        mod, "IVolumeSelector");

    using NameSelector = Acts::Geant4PhysicalVolumeSelectors::NameSelector;
    auto ns = py::class_<NameSelector, std::shared_ptr<NameSelector>>(
                  mod, "VolumeNameSelector", is)
                  .def(py::init<const std::vector<std::string>&, bool>());

    using Factory = Acts::Geant4DetectorSurfaceFactory;
    auto o = py::class_<Factory::Options>(mod, "SurfaceFactoryOptions")
                 .def(py::init<>());
    ACTS_PYTHON_STRUCT_BEGIN(o, Factory::Options);
    ACTS_PYTHON_MEMBER(scaleConversion);
    ACTS_PYTHON_MEMBER(convertMaterial);
    ACTS_PYTHON_MEMBER(convertedMaterialThickness);
    ACTS_PYTHON_MEMBER(sensitiveSurfaceSelector);
    ACTS_PYTHON_MEMBER(passiveSurfaceSelector);
    ACTS_PYTHON_STRUCT_END();
  }
  {
    py::class_<Acts::Geant4DetectorElement,
               std::shared_ptr<Acts::Geant4DetectorElement>>(
        mod, "Geant4DetectorElement");

    using Geant4Detector = ActsExamples::Geant4::Geant4Detector;

    auto g =
        py::class_<Geant4Detector, std::shared_ptr<Geant4Detector>>(
            mod, "Geant4Detector")
            .def(py::init<>())
            .def("constructDetector",
                 py::overload_cast<
                     const ActsExamples::Geant4::Geant4Detector::Config&>(
                     &ActsExamples::Geant4::Geant4Detector::constructDetector))
            .def("constructTrackingGeometry",
                 py::overload_cast<
                     const ActsExamples::Geant4::Geant4Detector::Config&>(
                     &ActsExamples::Geant4::Geant4Detector::
                         constructTrackingGeometry));

    auto c = py::class_<Geant4Detector::Config>(g, "Config").def(py::init<>());
    ACTS_PYTHON_STRUCT_BEGIN(c, Geant4Detector::Config);
    ACTS_PYTHON_MEMBER(name);
    ACTS_PYTHON_MEMBER(logLevel);
    ACTS_PYTHON_MEMBER(g4World);
    ACTS_PYTHON_MEMBER(g4SurfaceOptions);
    ACTS_PYTHON_MEMBER(protoDetector);
    ACTS_PYTHON_STRUCT_END();
  }

  Acts::Python::Context ctx;
  ctx.modules["geant4"] = &mod;

  addGeant4HepMC3(ctx);
}

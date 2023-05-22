// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Detector/Detector.hpp"
#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/MagneticField/MagneticFieldProvider.hpp"
#include "Acts/Plugins/Geant4/Geant4DetectorElement.hpp"
#include "Acts/Plugins/Geant4/Geant4DetectorSurfaceFactory.hpp"
#include "Acts/Plugins/Geant4/Geant4PhysicalVolumeSelectors.hpp"
#include "Acts/Plugins/Python/Utilities.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/Framework/IContextDecorator.hpp"
#include "ActsExamples/Geant4/ActsSteppingActionList.hpp"
#include "ActsExamples/Geant4/GdmlDetectorConstruction.hpp"
#include "ActsExamples/Geant4/Geant4Simulation.hpp"
#include "ActsExamples/Geant4/MagneticFieldWrapper.hpp"
#include "ActsExamples/Geant4/MaterialPhysicsList.hpp"
#include "ActsExamples/Geant4/MaterialSteppingAction.hpp"
#include "ActsExamples/Geant4/ParticleKillAction.hpp"
#include "ActsExamples/Geant4/ParticleTrackingAction.hpp"
#include "ActsExamples/Geant4/SensitiveSteppingAction.hpp"
#include "ActsExamples/Geant4/SensitiveSurfaceMapper.hpp"
#include "ActsExamples/Geant4/SimParticleTranslation.hpp"
#include "ActsExamples/Geant4Detector/Geant4Detector.hpp"
#include "ActsExamples/MuonSpectrometerMockupDetector/MockupSectorBuilder.hpp"
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
      randomNumbers, runManager, primaryGeneratorAction, runAction, eventAction,
      trackingAction, steppingAction, detectorConstruction, magneticField,
      sensitiveSurfaceMapper);

  auto makeGeant4Config =
      [](const Acts::Logger& logger,
         std::shared_ptr<const ActsExamples::RandomNumbers> randomNumbers,
         G4VUserDetectorConstruction* detector, G4VUserPhysicsList* physicsList,
         const SimParticleTranslation::Config& prCfg)
      -> Geant4Simulation::Config {
    Geant4Simulation::Config g4Cfg;

    // Set the main Geant4 algorithm, primary generation, detector
    // construction
    g4Cfg.randomNumbers = std::move(randomNumbers);
    g4Cfg.runManager = std::make_shared<G4RunManager>();
    g4Cfg.runManager->SetUserInitialization(physicsList);

    // Set the primarty generator
    g4Cfg.primaryGeneratorAction = new SimParticleTranslation(
        prCfg, logger.cloneWithSuffix("SimParticleTranslation"));
    g4Cfg.detectorConstruction = detector;

    return g4Cfg;
  };

  mod.def(
      "makeGeant4MaterialRecordingConfig",
      [makeGeant4Config](
          Acts::Logging::Level level, G4VUserDetectorConstruction* detector,
          std::shared_ptr<const ActsExamples::RandomNumbers> randomNumbers,
          const std::string& inputParticles,
          const std::string& outputMaterialTracks) {
        auto logger = Acts::getDefaultLogger("Geant4", level);
        auto physicsList = new MaterialPhysicsList(
            logger->cloneWithSuffix("MaterialPhysicsList"));

        // Read the particle from the generator
        SimParticleTranslation::Config g4PrCfg;
        g4PrCfg.forcedPdgCode = 0;
        g4PrCfg.forcedCharge = 0.;
        g4PrCfg.forcedMass = 0.;

        auto g4Cfg = makeGeant4Config(*logger, std::move(randomNumbers),
                                      detector, physicsList, g4PrCfg);
        g4Cfg.inputParticles = inputParticles;

        MaterialSteppingAction::Config mStepCfg;
        mStepCfg.excludeMaterials = {"Air", "Vacuum"};
        auto steppingAction = new MaterialSteppingAction(
            mStepCfg, logger->cloneWithSuffix("MaterialSteppingAction"));
        g4Cfg.steppingAction = steppingAction;

        // Set the material tracks at output
        g4Cfg.outputMaterialTracks = outputMaterialTracks;

        return g4Cfg;
      },
      "level"_a, "detector"_a, "randomNumbers"_a, "inputParticles"_a,
      "outputMaterialTracks"_a);

  mod.def(
      "makeGeant4SimulationConfig",
      [makeGeant4Config](
          Acts::Logging::Level level, G4VUserDetectorConstruction* detector,
          std::shared_ptr<const ActsExamples::RandomNumbers> randomNumbers,
          const std::string& inputParticles,
          const std::shared_ptr<const Acts::TrackingGeometry>& trackingGeometry,
          const std::shared_ptr<const Acts::MagneticFieldProvider>&
              magneticField,
          const std::vector<std::string>& volumeMappings,
          const std::vector<std::string>& materialMappings,
          std::shared_ptr<const Acts::Volume> killVolume,
          bool recordHitsOfSecondaries, bool keepParticlesWithoutHits) {
        auto logger = Acts::getDefaultLogger("Geant4", level);

        auto physicsList = new FTFP_BERT();
        auto g4Cfg =
            makeGeant4Config(*logger, std::move(randomNumbers), detector,
                             physicsList, SimParticleTranslation::Config{});
        g4Cfg.inputParticles = inputParticles;

        // Particle action
        ParticleTrackingAction::Config trackingCfg;
        trackingCfg.keepParticlesWithoutHits = keepParticlesWithoutHits;
        g4Cfg.trackingAction = new ParticleTrackingAction(
            trackingCfg, logger->cloneWithSuffix("ParticleTracking"));

        // Stepping actions
        ActsSteppingActionList::Config steppingCfg;

        SensitiveSteppingAction::Config g4StepCfg;
        g4StepCfg.charged = true;
        g4StepCfg.neutral = false;
        g4StepCfg.primary = true;
        g4StepCfg.secondary = recordHitsOfSecondaries;
        steppingCfg.actions.push_back(new SensitiveSteppingAction(
            g4StepCfg, logger->cloneWithSuffix("SensitiveStepping")));

        steppingCfg.actions.push_back(
            new ParticleKillAction(ParticleKillAction::Config{killVolume},
                                   logger->cloneWithSuffix("Killer")));

        g4Cfg.steppingAction = new ActsSteppingActionList(steppingCfg);

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
                  ssmCfg, logger->cloneWithSuffix("SensitiveSurfaceMapper"));
        }

        return g4Cfg;
      },
      "level"_a, "detector"_a, "randomNumbers"_a, "inputParticles"_a,
      py::arg("trackingGeometry") = nullptr, py::arg("magneticField") = nullptr,
      py::arg("volumeMappings") = std::vector<std::string>{},
      py::arg("materialMappings") = std::vector<std::string>{},
      py::arg("killVolume") = nullptr,
      py::arg("recordHitsOfSecondaries") = true,
      py::arg("keepParticlesWithoutHits") = true);

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
            .def(
                "constructDetector",
                [](Geant4Detector& self, const Geant4Detector::Config& cfg,
                   Logging::Level logLevel) {
                  auto logger = getDefaultLogger("Geant4Detector", logLevel);
                  return self.constructDetector(cfg, *logger);
                },
                py::arg("cfg"), py::arg("logLevel") = Logging::INFO)
            .def(
                "constructTrackingGeometry",
                [](Geant4Detector& self, const Geant4Detector::Config& cfg,
                   Logging::Level logLevel) {
                  auto logger = getDefaultLogger("Geant4Detector", logLevel);
                  return self.constructTrackingGeometry(cfg, *logger);
                },
                py::arg("cfg"), py::arg("logLevel") = Logging::INFO);

    auto c = py::class_<Geant4Detector::Config>(g, "Config").def(py::init<>());
    ACTS_PYTHON_STRUCT_BEGIN(c, Geant4Detector::Config);
    ACTS_PYTHON_MEMBER(name);
    ACTS_PYTHON_MEMBER(g4World);
    ACTS_PYTHON_MEMBER(g4SurfaceOptions);
    ACTS_PYTHON_MEMBER(protoDetector);
    ACTS_PYTHON_MEMBER(geometryIdentifierHook);
    ACTS_PYTHON_MEMBER(logLevel);
    ACTS_PYTHON_STRUCT_END();
  }

  {
    using MockupSectorBuilder = ActsExamples::MockupSectorBuilder;
    using Config = ActsExamples::MockupSectorBuilder::Config;
    using ChamberConfig = ActsExamples::MockupSectorBuilder::ChamberConfig;

    auto ms =
        py::class_<MockupSectorBuilder, std::shared_ptr<MockupSectorBuilder>>(
            mod, "MockupSectorBuilder")
            .def(py::init<const Config&>())
            .def("buildChamber", &MockupSectorBuilder::buildChamber)
            .def("buildSector", &MockupSectorBuilder::buildSector)
            .def("drawSector", &MockupSectorBuilder::drawSector);

    auto c = py::class_<Config>(ms, "Config").def(py::init<>());
    ACTS_PYTHON_STRUCT_BEGIN(c, Config);
    ACTS_PYTHON_MEMBER(gdmlPath);
    ACTS_PYTHON_MEMBER(NumberOfSectors);
    ACTS_PYTHON_MEMBER(toleranceOverlap);
    ACTS_PYTHON_STRUCT_END();

    auto cch = py::class_<ChamberConfig>(ms, "ChamberConfig").def(py::init<>());
    ACTS_PYTHON_STRUCT_BEGIN(cch, ChamberConfig);
    ACTS_PYTHON_MEMBER(name);
    ACTS_PYTHON_MEMBER(SensitiveNames);
    ACTS_PYTHON_MEMBER(PassiveNames);
    ACTS_PYTHON_STRUCT_END();
  }

  Acts::Python::Context ctx;
  ctx.modules["geant4"] = mod;

  addGeant4HepMC3(ctx);
}

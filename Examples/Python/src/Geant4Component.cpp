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

#include <G4RunManager.hh>
#include <G4UserEventAction.hh>
#include <G4UserRunAction.hh>
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
  py::class_<G4VUserDetectorConstruction,
             std::shared_ptr<G4VUserDetectorConstruction>>(
      mod, "G4VUserDetectorConstruction");

  py::class_<G4VPhysicalVolume, std::shared_ptr<G4VPhysicalVolume>>(
      mod, "G4VPhysicalVolume");

  // This is the actual class we're binding
  py::class_<GdmlDetectorConstruction, G4VUserDetectorConstruction,
             std::shared_ptr<GdmlDetectorConstruction>>(
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

  mod.def("makeGeant4MaterialRecordingConfig",
          Geant4Simulation::makeGeant4MaterialRecordingConfig, "level"_a,
          "detector"_a, "randomNumbers"_a, "inputParticles"_a,
          "outputMaterialTracks"_a);

  mod.def("makeGeant4SimulationConfig",
          Geant4Simulation::makeGeant4SimulationConfig, "level"_a, "detector"_a,
          "randomNumbers"_a, "inputParticles"_a,
          py::arg("trackingGeometry") = nullptr,
          py::arg("magneticField") = nullptr,
          py::arg("volumeMappings") = std::vector<std::string>{},
          py::arg("materialMappings") = std::vector<std::string>{},
          py::arg("killVolume") = nullptr,
          py::arg("killAfterTime") = std::numeric_limits<double>::infinity(),
          py::arg("recordHitsOfSecondaries") = true,
          py::arg("keepParticlesWithoutHits") = true);

  {
    using Detector = ActsExamples::Telescope::TelescopeDetector;
    using DetectorConstruction =
        ActsExamples::Telescope::TelescopeG4DetectorConstruction;

    py::class_<DetectorConstruction, G4VUserDetectorConstruction,
               std::shared_ptr<DetectorConstruction>>(
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

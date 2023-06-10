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
#include "ActsExamples/Geant4/DetectorConstructionFactory.hpp"
#include "ActsExamples/Geant4/GdmlDetectorConstruction.hpp"
#include "ActsExamples/Geant4/Geant4Simulation.hpp"
#include "ActsExamples/Geant4Detector/Geant4Detector.hpp"
#include "ActsExamples/MuonSpectrometerMockupDetector/MockupSectorBuilder.hpp"
#include "ActsExamples/TelescopeDetector/TelescopeG4DetectorConstruction.hpp"

#include <memory>

#include <G4RunManager.hh>
#include <G4UserEventAction.hh>
#include <G4UserRunAction.hh>
#include <pybind11/functional.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;
using namespace pybind11::literals;

using namespace ActsExamples;
using namespace Acts;
using namespace Acts::Python;

namespace Acts::Python {
void addGeant4HepMC3(Context& ctx);
}

PYBIND11_MODULE(ActsPythonBindingsGeant4, mod) {
  py::class_<DetectorConstructionFactory,
             std::shared_ptr<DetectorConstructionFactory>>(
      mod, "DetectorConstructionFactory");

  {
    using Algorithm = Geant4Simulation;
    using Config = Algorithm::Config;
    auto alg = py::class_<Algorithm, ActsExamples::IAlgorithm,
                          std::shared_ptr<Algorithm>>(mod, "Geant4Simulation")
                   .def(py::init<const Config&, Acts::Logging::Level>(),
                        py::arg("config"), py::arg("level"))
                   .def_property_readonly("config", &Algorithm::config);

    auto c1 = py::class_<Config, std::shared_ptr<Config>>(alg, "Config")
                  .def(py::init<>());
    ACTS_PYTHON_STRUCT_BEGIN(c1, Config);
    ACTS_PYTHON_MEMBER(inputParticles);
    ACTS_PYTHON_MEMBER(randomNumbers);
    ACTS_PYTHON_MEMBER(detectorConstructionFactory);
    ACTS_PYTHON_MEMBER(outputSimHits);
    ACTS_PYTHON_MEMBER(outputParticlesInitial);
    ACTS_PYTHON_MEMBER(outputParticlesFinal);
    ACTS_PYTHON_MEMBER(trackingGeometry);
    ACTS_PYTHON_MEMBER(magneticField);
    ACTS_PYTHON_MEMBER(volumeMappings);
    ACTS_PYTHON_MEMBER(materialMappings);
    ACTS_PYTHON_MEMBER(killVolume);
    ACTS_PYTHON_MEMBER(killAfterTime);
    ACTS_PYTHON_MEMBER(recordHitsOfSecondaries);
    ACTS_PYTHON_MEMBER(keepParticlesWithoutHits);
    ACTS_PYTHON_STRUCT_END();
  }

  {
    using Algorithm = Geant4MaterialRecording;
    using Config = Algorithm::Config;
    auto alg =
        py::class_<Algorithm, ActsExamples::IAlgorithm,
                   std::shared_ptr<Algorithm>>(mod, "Geant4MaterialRecording")
            .def(py::init<const Config&, Acts::Logging::Level>(),
                 py::arg("config"), py::arg("level"))
            .def_property_readonly("config", &Algorithm::config);

    auto c = py::class_<Config, std::shared_ptr<Config>>(alg, "Config")
                 .def(py::init<>());
    ACTS_PYTHON_STRUCT_BEGIN(c, Config);
    ACTS_PYTHON_MEMBER(inputParticles);
    ACTS_PYTHON_MEMBER(randomNumbers);
    ACTS_PYTHON_MEMBER(detectorConstructionFactory);
    ACTS_PYTHON_MEMBER(outputMaterialTracks);
    ACTS_PYTHON_STRUCT_END();
  }

  mod.def(
      "makeGeant4SimulationConfig",
      [](std::shared_ptr<DetectorConstructionFactory>
             detectorConstructionFactory,
         std::shared_ptr<const ActsExamples::RandomNumbers> randomNumbers,
         const std::string& inputParticles, const std::string& outputSimHits,
         const std::string& outputParticlesInitial,
         const std::string& outputParticlesFinal,
         std::shared_ptr<const Acts::TrackingGeometry> trackingGeometry,
         std::shared_ptr<const Acts::MagneticFieldProvider> magneticField,
         const std::vector<std::string>& volumeMappings,
         const std::vector<std::string>& materialMappings,
         std::shared_ptr<const Acts::Volume> killVolume, double killAfterTime,
         bool recordHitsOfSecondaries, bool keepParticlesWithoutHits) {
        Geant4Simulation::Config config;
        config.randomNumbers = std::move(randomNumbers);
        config.detectorConstructionFactory =
            std::move(detectorConstructionFactory);
        config.inputParticles = inputParticles;
        config.outputSimHits = outputSimHits;
        config.outputParticlesInitial = outputParticlesInitial;
        config.outputParticlesFinal = outputParticlesFinal;
        config.trackingGeometry = std::move(trackingGeometry);
        config.magneticField = std::move(magneticField);
        config.volumeMappings = volumeMappings;
        config.materialMappings = materialMappings;
        config.killVolume = std::move(killVolume);
        config.killAfterTime = killAfterTime;
        config.recordHitsOfSecondaries = recordHitsOfSecondaries;
        config.keepParticlesWithoutHits = keepParticlesWithoutHits;
        return config;
      },
      "detectorConstructionFactory"_a, "randomNumbers"_a, "inputParticles"_a,
      "outputSimHits"_a, "outputParticlesInitial"_a, "outputParticlesFinal"_a,
      py::arg("trackingGeometry") = nullptr, py::arg("magneticField") = nullptr,
      py::arg("volumeMappings") = std::vector<std::string>{},
      py::arg("materialMappings") = std::vector<std::string>{},
      py::arg("killVolume") = nullptr,
      py::arg("killAfterTime") = std::numeric_limits<double>::infinity(),
      py::arg("recordHitsOfSecondaries") = true,
      py::arg("keepParticlesWithoutHits") = true);

  mod.def(
      "makeGeant4MaterialRecordingConfig",
      [](std::shared_ptr<DetectorConstructionFactory>
             detectorConstructionFactory,
         std::shared_ptr<const ActsExamples::RandomNumbers> randomNumbers,
         const std::string& inputParticles,
         const std::string& outputMaterialTracks) {
        Geant4MaterialRecording::Config config;
        config.randomNumbers = std::move(randomNumbers);
        config.detectorConstructionFactory =
            std::move(detectorConstructionFactory);
        config.inputParticles = inputParticles;
        config.outputMaterialTracks = outputMaterialTracks;
        return config;
      },
      "detectorConstructionFactory"_a, "randomNumbers"_a, "inputParticles"_a,
      "outputMaterialTracks"_a);

  {
    py::class_<GdmlDetectorConstructionFactory, DetectorConstructionFactory,
               std::shared_ptr<GdmlDetectorConstructionFactory>>(
        mod, "GdmlDetectorConstructionFactory")
        .def(py::init<std::string>());
  }

  {
    py::class_<
        Telescope::TelescopeG4DetectorConstructionFactory,
        DetectorConstructionFactory,
        std::shared_ptr<Telescope::TelescopeG4DetectorConstructionFactory>>(
        mod, "TelescopeG4DetectorConstructionFactory")
        .def(py::init<const Telescope::TelescopeDetector::Config&>());
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

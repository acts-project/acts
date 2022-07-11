// This file is part of the Acts project.
//
// Copyright (C) 2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Plugins/Python/Utilities.hpp"
#include "ActsExamples/Io/EDM4hep/EDM4hepMeasurementReader.hpp"
#include "ActsExamples/Io/EDM4hep/EDM4hepMeasurementWriter.hpp"
#include "ActsExamples/Io/EDM4hep/EDM4hepMultiTrajectoryWriter.hpp"
#include "ActsExamples/Io/EDM4hep/EDM4hepParticleReader.hpp"
#include "ActsExamples/Io/EDM4hep/EDM4hepParticleWriter.hpp"
#include "ActsExamples/Io/EDM4hep/EDM4hepSimHitReader.hpp"
#include "ActsExamples/Io/EDM4hep/EDM4hepSimHitWriter.hpp"

#include <memory>

#include <pybind11/pybind11.h>

namespace py = pybind11;
using namespace pybind11::literals;

using namespace Acts;

namespace Acts::Python {

void addEDM4hep(Context& ctx) {
  auto mex = ctx.get("examples");

  auto edm4hep = mex.def_submodule("_edm4hep");

  {
    using Reader = ActsExamples::EDM4hepSimHitReader;
    using Config = Reader::Config;
    auto r = py::class_<Reader, ActsExamples::IReader, std::shared_ptr<Reader>>(
                 edm4hep, "EDM4hepSimHitReader")
                 .def(py::init<const Config&, Acts::Logging::Level>(),
                      py::arg("config"), py::arg("level"))
                 .def_property_readonly("config", &Reader::config);

    auto c = py::class_<Config>(r, "Config").def(py::init<>());
    ACTS_PYTHON_STRUCT_BEGIN(c, Config);
    ACTS_PYTHON_MEMBER(inputPath);
    ACTS_PYTHON_MEMBER(inputParticles);
    ACTS_PYTHON_MEMBER(outputSimHits);
    ACTS_PYTHON_MEMBER(dd4hepGeometryService);
    ACTS_PYTHON_STRUCT_END();
  }

  {
    using Writer = ActsExamples::EDM4hepSimHitWriter;
    using Config = Writer::Config;
    auto w = py::class_<Writer, ActsExamples::IWriter, std::shared_ptr<Writer>>(
                 edm4hep, "EDM4hepSimHitWriter")
                 .def(py::init<const Config&, Acts::Logging::Level>(),
                      py::arg("config"), py::arg("level"));

    auto c = py::class_<Config>(w, "Config").def(py::init<>());
    ACTS_PYTHON_STRUCT_BEGIN(c, Config);
    ACTS_PYTHON_MEMBER(inputSimHits);
    ACTS_PYTHON_MEMBER(inputParticles);
    ACTS_PYTHON_MEMBER(outputPath);
    ACTS_PYTHON_MEMBER(outputParticles);
    ACTS_PYTHON_MEMBER(outputSimTrackerHits);
    ACTS_PYTHON_STRUCT_END();
  }

  {
    using Reader = ActsExamples::EDM4hepMeasurementReader;
    using Config = Reader::Config;
    auto r = py::class_<Reader, ActsExamples::IReader, std::shared_ptr<Reader>>(
                 edm4hep, "EDM4hepMeasurementReader")
                 .def(py::init<const Config&, Acts::Logging::Level>(),
                      py::arg("config"), py::arg("level"))
                 .def_property_readonly("config", &Reader::config);

    auto c = py::class_<Config>(r, "Config").def(py::init<>());
    ACTS_PYTHON_STRUCT_BEGIN(c, Config);
    ACTS_PYTHON_MEMBER(inputPath);
    ACTS_PYTHON_MEMBER(outputMeasurements);
    ACTS_PYTHON_MEMBER(outputMeasurementSimHitsMap);
    ACTS_PYTHON_MEMBER(outputSourceLinks);
    ACTS_PYTHON_MEMBER(outputClusters);
    ACTS_PYTHON_STRUCT_END();
  }

  {
    using Writer = ActsExamples::EDM4hepMeasurementWriter;
    using Config = Writer::Config;
    auto w = py::class_<Writer, ActsExamples::IWriter, std::shared_ptr<Writer>>(
                 edm4hep, "EDM4hepMeasurementWriter")
                 .def(py::init<const Config&, Acts::Logging::Level>(),
                      py::arg("config"), py::arg("level"));

    auto c = py::class_<Config>(w, "Config").def(py::init<>());
    ACTS_PYTHON_STRUCT_BEGIN(c, Config);
    ACTS_PYTHON_MEMBER(inputMeasurements);
    ACTS_PYTHON_MEMBER(inputClusters);
    ACTS_PYTHON_MEMBER(inputSimHits);
    ACTS_PYTHON_MEMBER(inputMeasurementSimHitsMap);
    ACTS_PYTHON_MEMBER(outputPath);
    ACTS_PYTHON_STRUCT_END();
  }

  {
    using Reader = ActsExamples::EDM4hepParticleReader;
    using Config = Reader::Config;
    auto r = py::class_<Reader, ActsExamples::IReader, std::shared_ptr<Reader>>(
                 edm4hep, "EDM4hepParticleReader")
                 .def(py::init<const Config&, Acts::Logging::Level>(),
                      py::arg("config"), py::arg("level"))
                 .def_property_readonly("config", &Reader::config);

    auto c = py::class_<Config>(r, "Config").def(py::init<>());
    ACTS_PYTHON_STRUCT_BEGIN(c, Config);
    ACTS_PYTHON_MEMBER(inputPath);
    ACTS_PYTHON_MEMBER(inputParticles);
    ACTS_PYTHON_MEMBER(outputParticles);
    ACTS_PYTHON_STRUCT_END();
  }

  {
    using Writer = ActsExamples::EDM4hepParticleWriter;
    using Config = Writer::Config;
    auto w = py::class_<Writer, ActsExamples::IWriter, std::shared_ptr<Writer>>(
                 edm4hep, "EDM4hepParticleWriter")
                 .def(py::init<const Config&, Acts::Logging::Level>(),
                      py::arg("config"), py::arg("level"));

    auto c = py::class_<Config>(w, "Config").def(py::init<>());
    ACTS_PYTHON_STRUCT_BEGIN(c, Config);
    ACTS_PYTHON_MEMBER(inputParticles);
    ACTS_PYTHON_MEMBER(outputPath);
    ACTS_PYTHON_MEMBER(outputParticles);
    ACTS_PYTHON_STRUCT_END();
  }

  {
    using Writer = ActsExamples::EDM4hepMultiTrajectoryWriter;
    using Config = Writer::Config;
    auto w = py::class_<Writer, ActsExamples::IWriter, std::shared_ptr<Writer>>(
                 edm4hep, "EDM4hepMultiTrajectoryWriter")
                 .def(py::init<const Config&, Acts::Logging::Level>(),
                      py::arg("config"), py::arg("level"));

    auto c = py::class_<Config>(w, "Config").def(py::init<>());
    ACTS_PYTHON_STRUCT_BEGIN(c, Config);
    ACTS_PYTHON_MEMBER(inputTrajectories);
    ACTS_PYTHON_MEMBER(inputMeasurementParticlesMap);
    ACTS_PYTHON_MEMBER(outputPath);
    ACTS_PYTHON_STRUCT_END();
  }
}

}  // namespace Acts::Python

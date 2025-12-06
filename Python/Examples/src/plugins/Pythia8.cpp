// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/Generators/Pythia8ProcessGenerator.hpp"
#include "ActsPython/Utilities/Helpers.hpp"
#include "ActsPython/Utilities/Macros.hpp"

#include <memory>
#include <string>
#include <unordered_map>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/stl/filesystem.h>

namespace py = pybind11;
using namespace ActsExamples;
using namespace ActsPython;

PYBIND11_MODULE(ActsExamplesPythonBindingsPythia8, p8) {
  using Gen = Pythia8Generator;
  auto gen = py::class_<Gen, ParticlesGenerator, std::shared_ptr<Gen>>(
                 p8, "Pythia8Generator")
                 .def(py::init<const Gen::Config&, Acts::Logging::Level>(),
                      py::arg("config"), py::arg("level"));

  py::class_<Gen::Config>(gen, "Config")
      .def(py::init<>())
      .def_readwrite("pdgBeam0", &Gen::Config::pdgBeam0)
      .def_readwrite("pdgBeam1", &Gen::Config::pdgBeam1)
      .def_readwrite("cmsEnergy", &Gen::Config::cmsEnergy)
      .def_readwrite("settings", &Gen::Config::settings)
      .def_readwrite("printShortEventListing",
                     &Gen::Config::printShortEventListing)
      .def_readwrite("printLongEventListing",
                     &Gen::Config::printLongEventListing)
      .def_readwrite("labelSecondaries", &Gen::Config::labelSecondaries)
      .def_readwrite("spatialVertexThreshold",
                     &Gen::Config::spatialVertexThreshold)
      .def_readwrite("initializationSeed", &Gen::Config::initializationSeed)
      .def_readwrite("writeHepMC3", &Gen::Config::writeHepMC3);

  patchClassesWithConfig(p8);
}

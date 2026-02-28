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

#include <Pythia8/Pythia.h>
#include <pybind11/functional.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/stl/filesystem.h>

namespace py = pybind11;
using namespace ActsExamples;
using namespace ActsPython;

PYBIND11_MODULE(ActsExamplesPythonBindingsPythia8, p8) {
  // Pythia8 generator bindings
  {
    // Particle definition from Pythia8 with basic accessors, mother and
    // daughter info
    auto pythia8Particle =
        py::class_<Pythia8::Particle>(p8, "Pythia8Particle")
            .def("id", [](const Pythia8::Particle& self) { return self.id(); })
            .def("status",
                 [](const Pythia8::Particle& self) { return self.status(); })
            .def("px", [](const Pythia8::Particle& self) { return self.px(); })
            .def("py", [](const Pythia8::Particle& self) { return self.py(); })
            .def("pz", [](const Pythia8::Particle& self) { return self.pz(); })
            .def("pt",
                 [](const Pythia8::Particle& self) {
                   return std::sqrt(self.px() * self.px() +
                                    self.py() * self.py());
                 })
            .def("e", [](const Pythia8::Particle& self) { return self.e(); })
            .def("m", [](const Pythia8::Particle& self) { return self.m(); })
            .def("eta",
                 [](const Pythia8::Particle& self) { return self.eta(); })
            .def("mother1",
                 [](const Pythia8::Particle& self) { return self.mother1(); })
            .def("mother2",
                 [](const Pythia8::Particle& self) { return self.mother2(); })
            .def("daughter1",
                 [](const Pythia8::Particle& self) { return self.daughter1(); })
            .def("daughter2",
                 [](const Pythia8::Particle& self) { return self.daughter2(); })
            .def(
                "motherList",
                [](const Pythia8::Particle& self) { return self.motherList(); })
            .def("daughterList", [](const Pythia8::Particle& self) {
              return self.daughterList();
            });

    auto pythia8Event =
        py::class_<Pythia8::Event>(p8, "Pythia8Event")
            .def("size", [](const Pythia8::Event& self) { return self.size(); })
            .def("particles", [](const Pythia8::Event& self) {
              std::vector<Pythia8::Particle> particles;
              particles.reserve(self.size());
              for (int i = 0; i < self.size(); ++i) {
                particles.push_back(self[i]);
              }
              return particles;
            });

    auto pythia8 =
        py::class_<Pythia8::Pythia, std::shared_ptr<Pythia8::Pythia>>(p8,
                                                                      "Pythia8")
            .def(py::init<>())
            .def("event", [](Pythia8::Pythia& self) -> Pythia8::Event& {
              return self.event;
            });
  }

  // ACTS integrated Pythia8 generator
  {
    using Gen = Pythia8Generator;
    auto gen = py::class_<Gen, ParticlesGenerator, std::shared_ptr<Gen>>(
                   p8, "Pythia8Generator")
                   .def(py::init<const Gen::Config&, Acts::Logging::Level>(),
                        py::arg("config"), py::arg("level"));

    auto config = py::class_<Gen::Config>(gen, "Config").def(py::init<>());

    ACTS_PYTHON_STRUCT(config, pdgBeam0, pdgBeam1, cmsEnergy, settings,
                       printShortEventListing, printLongEventListing,
                       labelSecondaries, spatialVertexThreshold,
                       initializationSeed, writeHepMC3, eventSelectors);
  }
}

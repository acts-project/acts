// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Definitions/PdgParticle.hpp"
#include "Acts/EventData/ParticleHypothesis.hpp"
#include "ActsPython/Utilities/Helpers.hpp"

#include <type_traits>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;

using namespace Acts;

namespace ActsPython {

void addEventData(Context& ctx) {
  auto [m, mex] = ctx.get("main", "examples");

  py::class_<ParticleHypothesis>(m, "ParticleHypothesis")
      .def(py::init([](PdgParticle absPdg, float mass, float absCharge) {
             return ParticleHypothesis(absPdg, mass, AnyCharge{absCharge});
           }),
           py::arg("pdg"), py::arg("mass"), py::arg("absCharge"))
      .def(py::init([](std::underlying_type_t<PdgParticle> absPdg, float mass,
                       float absCharge) {
             return ParticleHypothesis(static_cast<PdgParticle>(absPdg), mass,
                                       AnyCharge{absCharge});
           }),
           py::arg("absPdg"), py::arg("mass"), py::arg("absCharge"))
      .def("__str__",
           [](const ParticleHypothesis& particleHypothesis) {
             std::stringstream os;
             particleHypothesis.toStream(os);
             return os.str();
           })
      .def("absolutePdg",
           [](const ParticleHypothesis& p) { return p.absolutePdg(); })
      .def("mass", [](const ParticleHypothesis& p) { return p.mass(); })
      .def("absoluteCharge",
           [](const ParticleHypothesis& p) { return p.absoluteCharge(); })
      .def_property_readonly_static(
          "muon",
          [](py::object /* self */) { return ParticleHypothesis::muon(); })
      .def_property_readonly_static(
          "pion",
          [](py::object /* self */) { return ParticleHypothesis::pion(); })
      .def_property_readonly_static(
          "electron",
          [](py::object /* self */) { return ParticleHypothesis::electron(); })
      .def_property_readonly_static(
          "kaon",
          [](py::object /* self */) { return ParticleHypothesis::kaon(); })
      .def_property_readonly_static(
          "proton",
          [](py::object /* self */) { return ParticleHypothesis::proton(); })
      .def_property_readonly_static(
          "geantino",
          [](py::object /* self */) { return ParticleHypothesis::geantino(); })
      .def_property_readonly_static(
          "chargedGeantino", [](py::object /* self */) {
            return ParticleHypothesis::chargedGeantino();
          });
}

}  // namespace ActsPython

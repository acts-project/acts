// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/EventData/ParticleHypothesis.hpp"
#include "Acts/Plugins/Python/Utilities.hpp"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;

using namespace Acts;

namespace Acts::Python {

void addEventData(Context& ctx) {
  auto [m, mex] = ctx.get("main", "examples");

  py::class_<Acts::ParticleHypothesis>(m, "ParticleHypothesis")
      .def("__str__",
           [](const Acts::ParticleHypothesis& particleHypothesis) {
             std::stringstream os;
             particleHypothesis.toStream(os);
             return os.str();
           })
      .def_property_readonly_static("muon",
                                    [](py::object /* self */) {
                                      return Acts::ParticleHypothesis::muon();
                                    })
      .def_property_readonly_static("pion",
                                    [](py::object /* self */) {
                                      return Acts::ParticleHypothesis::pion();
                                    })
      .def_property_readonly_static(
          "electron",
          [](py::object /* self */) {
            return Acts::ParticleHypothesis::electron();
          })
      .def_property_readonly_static(
          "geantino",
          [](py::object /* self */) {
            return Acts::ParticleHypothesis::geantino();
          })
      .def_property_readonly_static(
          "chargedGeantino", [](py::object /* self */) {
            return Acts::ParticleHypothesis::chargedGeantino();
          });
}

}  // namespace Acts::Python

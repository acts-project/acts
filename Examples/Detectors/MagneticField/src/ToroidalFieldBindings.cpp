// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/MagneticField/MagneticFieldProvider.hpp"
#include "ActsExamples/MagneticField/ToroidalField.hpp"

#include <memory>
#include <sstream>
#include <stdexcept>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;
using namespace Acts;

/// @brief Get the value of a field, throwing an exception if the result is invalid.
Vector3 getField(ToroidalField& self, const Vector3& position,
                 MagneticFieldProvider::Cache& cache) {
  if (Result<Vector3> res = self.getField(position, cache); !res.ok()) {
    std::stringstream ss;
    ss << "Field lookup failure with error: \"" << res.error() << "\"";
    throw std::runtime_error{ss.str()};
  } else {
    return *res;
  }
}

PYBIND11_MODULE(acts_toroidal_field, m) {
  m.doc() = "Toroidal Magnetic Field Python bindings";

  // BarrelConfig
  py::class_<ToroidalField::BarrelConfig>(m, "BarrelConfig")
      .def(py::init<>())
      .def_readwrite("R_in", &ToroidalField::BarrelConfig::R_in)
      .def_readwrite("R_out", &ToroidalField::BarrelConfig::R_out)
      .def_readwrite("c", &ToroidalField::BarrelConfig::c)
      .def_readwrite("b", &ToroidalField::BarrelConfig::b)
      .def_readwrite("I", &ToroidalField::BarrelConfig::I)
      .def_readwrite("Nturns", &ToroidalField::BarrelConfig::Nturns);

  // ECTConfig
  py::class_<ToroidalField::ECTConfig>(m, "ECTConfig")
      .def(py::init<>())
      .def_readwrite("R_in", &ToroidalField::ECTConfig::R_in)
      .def_readwrite("R_out", &ToroidalField::ECTConfig::R_out)
      .def_readwrite("c", &ToroidalField::ECTConfig::c)
      .def_readwrite("b", &ToroidalField::ECTConfig::b)
      .def_readwrite("I", &ToroidalField::ECTConfig::I)
      .def_readwrite("Nturns", &ToroidalField::ECTConfig::Nturns)
      .def_readwrite("gap", &ToroidalField::ECTConfig::gap);

  // LayoutConfig
  py::class_<ToroidalField::LayoutConfig>(m, "LayoutConfig")
      .def(py::init<>())
      .def_readwrite("theta0_deg", &ToroidalField::LayoutConfig::theta0_deg)
      .def_readwrite("thetaStep_deg",
                     &ToroidalField::LayoutConfig::thetaStep_deg)
      .def_readwrite("nCoils", &ToroidalField::LayoutConfig::nCoils)
      .def_readwrite("nArc", &ToroidalField::LayoutConfig::nArc)
      .def_readwrite("nStraight", &ToroidalField::LayoutConfig::nStraight)
      .def_readwrite("closeLoop", &ToroidalField::LayoutConfig::closeLoop)
      .def_readwrite("eps", &ToroidalField::LayoutConfig::eps);

  // Main Config
  py::class_<ToroidalField::Config>(m, "Config")
      .def(py::init<>())
      .def_readwrite("barrel", &ToroidalField::Config::barrel)
      .def_readwrite("ect", &ToroidalField::Config::ect)
      .def_readwrite("layout", &ToroidalField::Config::layout)
      .def_readwrite("barrelSigns", &ToroidalField::Config::barrelSigns)
      .def_readwrite("ectSigns", &ToroidalField::Config::ectSigns);

  // ToroidalField class
  py::class_<ToroidalField, std::shared_ptr<ToroidalField>,
             Acts::MagneticFieldProvider>(m, "ToroidalField")
      .def(py::init<const ToroidalField::Config&>(), py::arg("config"))
      .def("getField", &getField)
      .def("makeCache", &ToroidalField::makeCache)
      .def("config", &ToroidalField::config,
           py::return_value_policy::reference_internal);
}

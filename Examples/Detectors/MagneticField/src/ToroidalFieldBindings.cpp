// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/MagneticField/MagneticFieldProvider.hpp"
#include "ActsExamples/MagneticField/ToroidalField.hpp"
#include "ActsExamples/MagneticField/ToroidalFieldMap.hpp"

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

  // EctConfig
  py::class_<ToroidalField::EctConfig>(m, "EctConfig")
      .def(py::init<>())
      .def_readwrite("R_in", &ToroidalField::EctConfig::R_in)
      .def_readwrite("R_out", &ToroidalField::EctConfig::R_out)
      .def_readwrite("c", &ToroidalField::EctConfig::c)
      .def_readwrite("b", &ToroidalField::EctConfig::b)
      .def_readwrite("I", &ToroidalField::EctConfig::I)
      .def_readwrite("Nturns", &ToroidalField::EctConfig::Nturns)
      .def_readwrite("gap", &ToroidalField::EctConfig::gap);

  // LayoutConfig
  py::class_<ToroidalField::LayoutConfig>(m, "LayoutConfig")
      .def(py::init<>())
      .def_readwrite("theta0", &ToroidalField::LayoutConfig::theta0)
      .def_readwrite("thetaStep", &ToroidalField::LayoutConfig::thetaStep)
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

  // ToroidalFieldMap functions
  m.def("toroidalFieldMapCyl", &Acts::toroidalFieldMapCyl, py::arg("rLim"),
        py::arg("phiLim"), py::arg("zLim"), py::arg("nBins"), py::arg("field"),
        "Create cylindrical toroidal field map");

  m.def("toroidalFieldMapXYZ", &Acts::toroidalFieldMapXYZ, py::arg("xLim"),
        py::arg("yLim"), py::arg("zLim"), py::arg("nBins"), py::arg("field"),
        "Create Cartesian toroidal field map");
}

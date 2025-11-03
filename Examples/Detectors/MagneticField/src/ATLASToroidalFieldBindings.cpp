// This file is part of the ACTS project.
//
// Copyright (C) 2025 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/MagneticField/MagneticFieldProvider.hpp"
#include "ActsExamples/MagneticField/ATLASToroidalField.hpp"

#include <memory>
#include <sstream>
#include <stdexcept>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;
using namespace Acts;

/// @brief Get the value of a field, throwing an exception if the result is invalid.
Vector3 getField(ATLASToroidalField& self, const Vector3& position,
                 MagneticFieldProvider::Cache& cache) {
  if (Result<Vector3> res = self.getField(position, cache); !res.ok()) {
    std::stringstream ss;
    ss << "Field lookup failure with error: \"" << res.error() << "\"";
    throw std::runtime_error{ss.str()};
  } else {
    return *res;
  }
}

PYBIND11_MODULE(acts_atlas_toroidal_field, m) {
  m.doc() = "ATLAS Toroidal Magnetic Field Python bindings";

  // BarrelConfig
  py::class_<ATLASToroidalField::BarrelConfig>(m, "BarrelConfig")
      .def(py::init<>())
      .def_readwrite("R_in", &ATLASToroidalField::BarrelConfig::R_in)
      .def_readwrite("R_out", &ATLASToroidalField::BarrelConfig::R_out)
      .def_readwrite("c", &ATLASToroidalField::BarrelConfig::c)
      .def_readwrite("b", &ATLASToroidalField::BarrelConfig::b)
      .def_readwrite("I", &ATLASToroidalField::BarrelConfig::I)
      .def_readwrite("Nturns", &ATLASToroidalField::BarrelConfig::Nturns);

  // ECTConfig
  py::class_<ATLASToroidalField::ECTConfig>(m, "ECTConfig")
      .def(py::init<>())
      .def_readwrite("R_in", &ATLASToroidalField::ECTConfig::R_in)
      .def_readwrite("R_out", &ATLASToroidalField::ECTConfig::R_out)
      .def_readwrite("c", &ATLASToroidalField::ECTConfig::c)
      .def_readwrite("b", &ATLASToroidalField::ECTConfig::b)
      .def_readwrite("I", &ATLASToroidalField::ECTConfig::I)
      .def_readwrite("Nturns", &ATLASToroidalField::ECTConfig::Nturns)
      .def_readwrite("gap", &ATLASToroidalField::ECTConfig::gap);

  // LayoutConfig
  py::class_<ATLASToroidalField::LayoutConfig>(m, "LayoutConfig")
      .def(py::init<>())
      .def_readwrite("theta0_deg",
                     &ATLASToroidalField::LayoutConfig::theta0_deg)
      .def_readwrite("thetaStep_deg",
                     &ATLASToroidalField::LayoutConfig::thetaStep_deg)
      .def_readwrite("nCoils", &ATLASToroidalField::LayoutConfig::nCoils)
      .def_readwrite("nArc", &ATLASToroidalField::LayoutConfig::nArc)
      .def_readwrite("nStraight", &ATLASToroidalField::LayoutConfig::nStraight)
      .def_readwrite("closeLoop", &ATLASToroidalField::LayoutConfig::closeLoop)
      .def_readwrite("eps", &ATLASToroidalField::LayoutConfig::eps);

  // Main Config
  py::class_<ATLASToroidalField::Config>(m, "Config")
      .def(py::init<>())
      .def_readwrite("barrel", &ATLASToroidalField::Config::barrel)
      .def_readwrite("ect", &ATLASToroidalField::Config::ect)
      .def_readwrite("layout", &ATLASToroidalField::Config::layout)
      .def_readwrite("barrelSigns", &ATLASToroidalField::Config::barrelSigns)
      .def_readwrite("ectSigns", &ATLASToroidalField::Config::ectSigns);

  // ATLASToroidalField class
  py::class_<ATLASToroidalField, std::shared_ptr<ATLASToroidalField>,
             Acts::MagneticFieldProvider>(m, "ATLASToroidalField")
      .def(py::init<const ATLASToroidalField::Config&>(), py::arg("config"))
      .def("getField", &getField)
      .def("makeCache", &ATLASToroidalField::makeCache)
      .def("config", &ATLASToroidalField::config,
           py::return_value_policy::reference_internal);
}

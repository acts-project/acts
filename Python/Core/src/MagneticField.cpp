// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/MagneticField/BFieldMapUtils.hpp"
#include "Acts/MagneticField/ConstantBField.hpp"
#include "Acts/MagneticField/MagneticFieldContext.hpp"
#include "Acts/MagneticField/MagneticFieldProvider.hpp"
#include "Acts/MagneticField/MultiRangeBField.hpp"
#include "Acts/MagneticField/NullBField.hpp"
#include "Acts/MagneticField/SolenoidBField.hpp"
#include "Acts/Utilities/Result.hpp"
#include "ActsPython/Utilities/Context.hpp"

#include <stdexcept>

#include <pybind11/eval.h>
#include <pybind11/pybind11.h>
#include <pybind11/pytypes.h>
#include <pybind11/stl.h>

namespace py = pybind11;
using namespace pybind11::literals;
using namespace Acts;

namespace ActsPython {

/// @brief Get the value of a field, throwing an exception if the result is
/// invalid.
Vector3 getField(MagneticFieldProvider& self, const Vector3& position,
                 MagneticFieldProvider::Cache& cache) {
  if (Result<Vector3> res = self.getField(position, cache); !res.ok()) {
    std::stringstream ss;

    ss << "Field lookup failure with error: \"" << res.error() << "\"";

    throw std::runtime_error{ss.str()};
  } else {
    return *res;
  }
}

/// This adds the definitions from Core/MagneticField to the python module
/// @param ctx the context container for the python modules
void addMagneticField(Context& ctx) {
  auto& m = ctx.get("main");

  // Add the MagneticFieldContext class to the module
  py::class_<MagneticFieldContext>(m, "MagneticFieldContext").def(py::init<>());

  // Add the MagneticFieldProvider and the fields
  py::class_<MagneticFieldProvider, std::shared_ptr<MagneticFieldProvider>>(
      m, "MagneticFieldProvider")
      .def("getField", &getField)
      .def("makeCache", &MagneticFieldProvider::makeCache);

  py::class_<InterpolatedMagneticField,
             std::shared_ptr<InterpolatedMagneticField>>(
      m, "InterpolatedMagneticField");

  m.def("solenoidFieldMap", &solenoidFieldMap, py::arg("rlim"), py::arg("zlim"),
        py::arg("nbins"), py::arg("field"));

  py::class_<ConstantBField, MagneticFieldProvider,
             std::shared_ptr<ConstantBField>>(m, "ConstantBField")
      .def(py::init<Vector3>());

  py::class_<NullBField, MagneticFieldProvider, std::shared_ptr<NullBField>>(
      m, "NullBField")
      .def(py::init<>());

  py::class_<MultiRangeBField, MagneticFieldProvider,
             std::shared_ptr<MultiRangeBField>>(m, "MultiRangeBField")
      .def(py::init<std::vector<std::pair<RangeXD<3, double>, Vector3>>>());

  {
    using Config = SolenoidBField::Config;

    auto sol = py::class_<SolenoidBField, MagneticFieldProvider,
                          std::shared_ptr<SolenoidBField>>(m, "SolenoidBField")
                   .def(py::init<Config>())
                   .def(py::init([](double radius, double length,
                                    std::size_t nCoils, double bMagCenter) {
                          return SolenoidBField{
                              Config{radius, length, nCoils, bMagCenter}};
                        }),
                        py::arg("radius"), py::arg("length"), py::arg("nCoils"),
                        py::arg("bMagCenter"));

    py::class_<Config>(sol, "Config")
        .def(py::init<>())
        .def_readwrite("radius", &Config::radius)
        .def_readwrite("length", &Config::length)
        .def_readwrite("nCoils", &Config::nCoils)
        .def_readwrite("bMagCenter", &Config::bMagCenter);
  }
}
}  // namespace ActsPython

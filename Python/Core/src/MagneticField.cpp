// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Definitions/Units.hpp"
#include "Acts/MagneticField/BFieldMapUtils.hpp"
#include "Acts/MagneticField/ConstantBField.hpp"
#include "Acts/MagneticField/MagneticFieldContext.hpp"
#include "Acts/MagneticField/MagneticFieldProvider.hpp"
#include "Acts/MagneticField/MultiRangeBField.hpp"
#include "Acts/MagneticField/NullBField.hpp"
#include "Acts/MagneticField/SolenoidBField.hpp"
#include "Acts/MagneticField/TextMagneticFieldIo.hpp"
#include "ActsPython/Utilities/Helpers.hpp"
#include "ActsPython/Utilities/Macros.hpp"

#include <array>
#include <cstddef>
#include <filesystem>
#include <memory>
#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;
using namespace pybind11::literals;

using namespace Acts;
using namespace ActsExamples;

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

/// @brief Add the magnetic field bindings to a module.
/// @param m the module to add the bindings to
void addMagneticField(py::module_& m) {
  py::class_<Acts::MagneticFieldContext>(m, "MagneticFieldContext")
      .def(py::init<>());

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

  using InterpolatedMagneticField2 = InterpolatedBFieldMap<
      Grid<Vector2, Axis<AxisType::Equidistant>, Axis<AxisType::Equidistant>>>;

  using InterpolatedMagneticField3 = InterpolatedBFieldMap<
      Grid<Vector3, Axis<AxisType::Equidistant>, Axis<AxisType::Equidistant>,
           Axis<AxisType::Equidistant>>>;

  py::class_<InterpolatedMagneticField2, InterpolatedMagneticField,
             MagneticFieldProvider,
             std::shared_ptr<InterpolatedMagneticField2>>(
      m, "InterpolatedMagneticField2");

  py::class_<InterpolatedMagneticField3, InterpolatedMagneticField,
             MagneticFieldProvider,
             std::shared_ptr<InterpolatedMagneticField3>>(
      m, "InterpolatedMagneticField3");

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

  m.def(
      "MagneticFieldMapXyz",
      [](const std::string& filename, double lengthUnit, double BFieldUnit,
         bool firstOctant) {
        const std::filesystem::path file = filename;

        auto mapBins = [](std::array<std::size_t, 3> bins,
                          std::array<std::size_t, 3> sizes) {
          return (bins[0] * (sizes[1] * sizes[2]) + bins[1] * sizes[2] +
                  bins[2]);
        };

        if (file.extension() == ".txt") {
          auto map = makeMagneticFieldMapXyzFromText(std::move(mapBins),
                                                     file.native(), lengthUnit,
                                                     BFieldUnit, firstOctant);
          return std::make_shared<decltype(map)>(std::move(map));
        } else {
          throw std::runtime_error("Unsupported magnetic field map file type");
        }
      },
      py::arg("file"), py::arg("lengthUnit") = UnitConstants::mm,
      py::arg("BFieldUnit") = UnitConstants::T, py::arg("firstOctant") = false);

  m.def(
      "MagneticFieldMapRz",
      [](const std::string& filename, double lengthUnit, double BFieldUnit,
         bool firstQuadrant) {
        const std::filesystem::path file = filename;

        auto mapBins = [](std::array<std::size_t, 2> bins,
                          std::array<std::size_t, 2> sizes) {
          return (bins[1] * sizes[0] + bins[0]);
        };

        if (file.extension() == ".txt") {
          auto map = makeMagneticFieldMapRzFromText(std::move(mapBins),
                                                    file.native(), lengthUnit,
                                                    BFieldUnit, firstQuadrant);
          return std::make_shared<decltype(map)>(std::move(map));
        } else {
          throw std::runtime_error("Unsupported magnetic field map file type");
        }
      },
      py::arg("file"), py::arg("lengthUnit") = UnitConstants::mm,
      py::arg("BFieldUnit") = UnitConstants::T,
      py::arg("firstQuadrant") = false);
}

}  // namespace ActsPython

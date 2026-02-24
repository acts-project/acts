// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsPlugins/Detray/DetrayConversionUtils.hpp"
#include "ActsPython/Utilities/Helpers.hpp"
#include "ActsPython/Utilities/Macros.hpp"

#include <memory>
#include <string>

#include <detray/builders/detector_builder.hpp>
#include <detray/geometry/barcode.hpp>
#include <detray/geometry/surface_descriptor.hpp>
#include <pybind11/operators.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <vecmem/memory/host_memory_resource.hpp>
#include <vecmem/memory/memory_resource.hpp>

namespace py = pybind11;
using namespace pybind11::literals;

PYBIND11_MODULE(ActsPluginsPythonBindingsDetray, detray) {
  using namespace ActsPlugins;

  {
    py::class_<DetrayHostDetector, std::shared_ptr<DetrayHostDetector>>(
        detray, "detray_detector");
  }

  using surface_id = detray::surface_id;
  py::enum_<surface_id>(detray, "surface_id")
      .value("e_portal", surface_id::e_portal)
      .value("e_sensitive", surface_id::e_sensitive)
      .value("e_passive", surface_id::e_passive)
      .value("e_unknown", surface_id::e_unknown)
      .value("e_all", surface_id::e_all);

  using barcode = detray::geometry::barcode;
  auto geometry = detray.def_submodule("geometry");
  py::class_<barcode>(geometry, "barcode")
      .def(py::init<>())
      .def(py::init<detray::geometry::barcode::value_t>())
      .def_property_readonly("value", &barcode::value)
      .def_property("volume", &barcode::volume, &barcode::set_volume)
      .def_property("id", &barcode::id, &barcode::set_id)
      .def_property("index", &barcode::index, &barcode::set_index)
      .def_property("transform", &barcode::transform, &barcode::set_transform)
      .def_property("extra", &barcode::extra, &barcode::set_extra)
      .def_property_readonly("is_invalid", &barcode::is_invalid)
      .def(py::self == py::self);

  // py::class_<detray::dindex>(detray, "dindex");

  using dtyped_index = detray::dtyped_index<detray::dindex, detray::dindex>;
  py::class_<detray::dtyped_index<detray::dindex, detray::dindex>>(
      detray, "dtyped_index")
      .def(py::init<>())
      .def_property("id", &dtyped_index::id, &dtyped_index::set_id)
      .def(py::init<const detray::dindex, const detray::dindex>())
      .def_property("index", &dtyped_index::index, &dtyped_index::set_index)
      .def_property("id", &dtyped_index::id, &dtyped_index::set_id)
      .def("shift", &dtyped_index::shift<dtyped_index::index_type>)
      .def_property_readonly("is_invalid", &dtyped_index::is_invalid)
      .def_property_readonly("is_invalid_id", &dtyped_index::is_invalid_id)
      .def_property_readonly("is_invalid_index",
                             &dtyped_index::is_invalid_index)
      .def(py::self == py::self)
      .def(py::self + py::self)
      .def(py::self += py::self)
      .def(py::self - py::self)
      .def(py::self -= py::self)
      .def("__str__", [](const dtyped_index& self) {
        std::stringstream ss;
        ss << self;
        return ss.str();
      });
  ;
}

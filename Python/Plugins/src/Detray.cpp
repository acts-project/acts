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
#include <detray/definitions/algebra.hpp>
#include <detray/detectors/toy_metadata.hpp>
#include <detray/geometry/barcode.hpp>
#include <detray/geometry/surface_descriptor.hpp>
#include <pybind11/operators.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <vecmem/memory/host_memory_resource.hpp>
#include <vecmem/memory/memory_resource.hpp>

namespace py = pybind11;
using namespace pybind11::literals;

template <typename metadata_t>
void addGeometry(py::module& m) {
  using metadata = metadata_t;
  using mask_link = metadata::mask_link;

  auto pyMetadata = py::class_<metadata>(m, "metadata");
  // auto pyMaskLink = py::class_<mask_link>(pyMetadata, "mask_link");

  py::enum_<typename metadata::mask_ids>(pyMetadata, "mask_ids")
      // py::enum_<typename mask_link::id_type>(pyMaskLink, "id_type")
      .value("e_rectangle2", metadata::mask_ids::e_rectangle2)
      .value("e_trapezoid2", metadata::mask_ids::e_trapezoid2)
      .value("e_portal_cylinder2", metadata::mask_ids::e_portal_cylinder2)
      .value("e_portal_ring2", metadata::mask_ids::e_portal_ring2)
      .value("e_cylinder2", metadata::mask_ids::e_cylinder2);

  auto pyMaskLink = py::class_<mask_link>(m, "mask_link")
                        .def(py::init<>())
                        .def(py::init<const typename mask_link::id_type,
                                      const typename mask_link::index_type>())
                        .def_property("id", &mask_link::id, &mask_link::set_id);

  py::class_<typename mask_link::index_type>(pyMaskLink, "index_type")
      .def(py::init<>());
}

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

  auto toy = detray.def_submodule("toy");
  addGeometry<detray::toy_metadata<detray::array<float>>>(toy);
}

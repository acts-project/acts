// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// Detray core include(s)
#include "detray/core/detector.hpp"

// Detray IO include(s)
#include "detray/io/frontend/detector_reader.hpp"
#include "detray/io/frontend/detector_reader_config.hpp"

// Detray algebra plugin + detector metadata
#include "algebra/array.hpp"
#include "detray/definitions/algebra.hpp"
#include "detray/detectors/default_metadata.hpp"

// Vecmem include(s)
#include <vecmem/memory/host_memory_resource.hpp>

// Pybind11 include(s)
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

// System include(s)
#include <memory>
#include <sstream>
#include <string>
#include <string_view>
#include <utility>

namespace py = pybind11;

#define STRINGIFY(x) #x
/// Expand macro and turn its value into a string.
#define STRINGIFY_HELPER(x) STRINGIFY(x)

namespace {

using algebra_t = detray::array<DETRAY_CUSTOM_SCALARTYPE>;
using detector_t = detray::detector<detray::default_metadata<algebra_t>>;
using volume_descriptor_t = detector_t::volume_type;
using volume_container_t = detector_t::volume_container;
using surface_descriptor_t = detector_t::surface_type;
using surface_store_t = detector_t::surface_lookup_container;
using surface_container_t = detector_t::surface_container;
using geometry_context_t = detector_t::geometry_context;
using transform_store_t = detector_t::transform_container;
using mask_store_t = detector_t::mask_container;
using material_store_t = detector_t::material_container;
using accelerator_store_t = detector_t::accelerator_container;

/// Owns a detector together with the memory resource its data lives in.
struct detector_handle {
  std::unique_ptr<vecmem::host_memory_resource> memory_resource;
  detector_t detector;
};

/// Read a detector (default metadata) from a JSON file.
std::pair<detector_handle, detray::name_map> read_detector(
    const std::string &file_name) {
  auto mr = std::make_unique<vecmem::host_memory_resource>();

  detray::io::detector_reader_config cfg{};
  cfg.add_file(file_name);

  auto [det, names] = detray::io::read_detector<detector_t>(*mr, cfg);

  return {detector_handle{std::move(mr), std::move(det)}, std::move(names)};
}

}  // namespace

PYBIND11_MODULE(DetrayPythonBindings, m) {
  m.doc() = "Detray core bindings";

  py::class_<volume_descriptor_t>(m, "VolumeDescriptor");
  py::class_<surface_descriptor_t>(m, "SurfaceDescriptor");
  py::class_<surface_store_t>(m, "SurfaceStore");
  py::class_<geometry_context_t>(m, "GeometryContext");
  py::class_<transform_store_t>(m, "TransformStore");
  py::class_<mask_store_t>(m, "MaskStore");
  py::class_<material_store_t>(m, "MaterialStore");
  py::class_<accelerator_store_t>(m, "AcceleratorStore");
  py::class_<detray::name_map>(m, "NameMap");

  py::class_<detector_handle>(
      m, "DetectorDefaultMetadata" STRINGIFY_HELPER(DETRAY_CUSTOM_SCALARTYPE))
      .def(
          "n_volumes",
          [](const detector_handle &d) { return d.detector.volumes().size(); },
          "Number of volumes in the detector")
      .def(
          "n_surfaces",
          [](const detector_handle &d) { return d.detector.surfaces().size(); },
          "Number of surfaces in the detector")
      .def(
          "name",
          [](const detector_handle &d, const detray::name_map &names) {
            return d.detector.name(names);
          },
          py::arg("names"), "Detector name")
      .def_property_readonly(
          "volumes",
          // Converts to list[VolumeDescriptor].
          [](const detector_handle &d) -> const volume_container_t & {
            return d.detector.volumes();
          },
          py::return_value_policy::reference_internal, "All volumes")
      .def_property_readonly(
          "surfaces",
          [](const detector_handle &d) -> const surface_store_t & {
            return d.detector.surfaces();
          },
          py::return_value_policy::reference_internal, "All surfaces")
      .def(
          "surface",
          [](const detector_handle &d,
             detray::dindex index) -> const surface_descriptor_t & {
            return d.detector.surface(index);
          },
          py::arg("index"), py::return_value_policy::reference_internal,
          "Surface by index")
      .def_property_readonly(
          "portals",
          // Converts to list[SurfaceDescriptor].
          [](const detector_handle &d) -> const surface_container_t & {
            return d.detector.portals();
          },
          "All portals")
      .def_property_readonly(
          "transformStore",
          [](const detector_handle &d) -> const transform_store_t & {
            return d.detector.transform_store();
          },
          py::return_value_policy::reference_internal, "Transform store")
      .def_property_readonly(
          "maskStore",
          [](const detector_handle &d) -> const mask_store_t & {
            return d.detector.mask_store();
          },
          py::return_value_policy::reference_internal, "Mask store")
      .def_property_readonly(
          "materialStore",
          [](const detector_handle &d) -> const material_store_t & {
            return d.detector.material_store();
          },
          py::return_value_policy::reference_internal, "Material store")
      .def_property_readonly(
          "acceleratorStore",
          [](const detector_handle &d) -> const accelerator_store_t & {
            return d.detector.accelerator_store();
          },
          py::return_value_policy::reference_internal, "Accelerator store")
      .def("__repr__", [](const detector_handle &d) {
        std::ostringstream os;
        os << d.detector;
        return os.str();
      });

  m.def("readDetector", &read_detector, py::arg("fileName"),
        "Read a detector from a JSON file");
}

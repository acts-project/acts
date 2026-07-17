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
#include <cstddef>
#include <memory>
#include <sstream>
#include <stdexcept>
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
  py::class_<surface_store_t>(m, "SurfaceStore")
      .def("__len__", [](const surface_store_t &s) { return s.size(); })
      .def(
          "empty", [](const surface_store_t &s) { return s.empty(); },
          "Whether the collection has no surfaces")
      .def(
          "__getitem__",
          [](const surface_store_t &s,
             py::ssize_t i) -> const surface_descriptor_t & {
            const auto n = static_cast<py::ssize_t>(s.size());
            if (i < 0) {
              i += n;
            }
            if (i < 0 || i >= n) {
              throw py::index_error();
            }
            return s[static_cast<std::size_t>(i)];
          },
          py::arg("index"), py::return_value_policy::reference_internal,
          "Surface at a given index")
      .def(
          "__iter__",
          [](const surface_store_t &s) {
            return py::make_iterator<
                py::return_value_policy::reference_internal,
                decltype(s.begin()), decltype(s.end()),
                const surface_descriptor_t &>(s.begin(), s.end());
          },
          py::keep_alive<0, 1>(), "Iterate over the surfaces");
  py::class_<geometry_context_t>(m, "GeometryContext");
  py::class_<transform_store_t>(m, "TransformStore");
  py::class_<mask_store_t>(m, "MaskStore");
  py::class_<material_store_t>(m, "MaterialStore");
  py::class_<accelerator_store_t>(m, "AcceleratorStore");
  py::class_<detray::name_map>(m, "NameMap")
      .def(py::init<>())
      .def_property(
          "detector_name",
          [](const detray::name_map &n) { return n.get_detector_name(); },
          [](detray::name_map &n, std::string_view name) {
            n.set_detector_name(name);
          },
          "Name of the detector")
      .def("empty", &detray::name_map::empty, "Whether no volume names are mapped")
      .def(
          "__contains__",
          [](const detray::name_map &n, detray::dindex index) {
            return n.contains(index);
          },
          py::arg("index"), "Whether a volume index is mapped")
      .def(
          "__contains__",
          [](const detray::name_map &n, std::string_view name) {
            return n.contains(name);
          },
          py::arg("name"), "Whether a volume name is mapped")
      .def(
          "__setitem__",
          [](detray::name_map &n, detray::dindex index, const std::string &name) {
            n.emplace(index, name);
          },
          py::arg("index"), py::arg("name"), "Map a volume index to a name")
      .def(
          "__getitem__",
          [](const detray::name_map &n, detray::dindex index) -> std::string {
            try {
              return n.at(index);
            } catch (const std::out_of_range &) {
              throw py::key_error(std::to_string(index));
            }
          },
          py::arg("index"), "Volume name at a volume index")
      .def(
          "__getitem__",
          [](const detray::name_map &n,
             std::string_view name) -> detray::dindex {
            try {
              return n.at(name);
            } catch (const std::out_of_range &) {
              throw py::key_error(std::string{name});
            }
          },
          py::arg("name"), "Volume index at a volume name")
      .def("clear", &detray::name_map::clear, "Clear detector and volume names")
      .def("clear_names", &detray::name_map::clear_names,
           "Clear volume names, keep the detector name")
      .def("__repr__", [](const detray::name_map &n) {
        return "NameMap(detector_name='" + n.get_detector_name() + "')";
      });

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

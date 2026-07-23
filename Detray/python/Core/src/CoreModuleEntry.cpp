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

// Detray propagation include(s)
#include "detray/propagator/propagation_config.hpp"

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

using reader_config_t = detray::io::detector_reader_config;
using intersection_config_t = detray::intersection::config;
using navigation_config_t = detray::navigation::config;
using stepping_config_t = detray::stepping::config;
using propagation_config_t = detray::propagation::config;

template <typename T>
std::string to_string(const T &obj) {
  std::ostringstream os;
  os << obj;
  return os.str();
}

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

/// Read a detector (default metadata) as configured by @param cfg .
std::pair<detector_handle, detray::name_map> read_detector(
    const reader_config_t &cfg) {
  auto mr = std::make_unique<vecmem::host_memory_resource>();

  auto [det, names] = detray::io::read_detector<detector_t>(*mr, cfg);

  return {detector_handle{std::move(mr), std::move(det)}, std::move(names)};
}

}  // namespace

PYBIND11_MODULE(DetrayPythonBindings, m) {
  m.doc() = "Detray core bindings";

  py::class_<reader_config_t>(m, "DetectorReaderConfig")
      .def(py::init<>())
      .def_property_readonly("files", &reader_config_t::files, "Input files")
      .def_property(
          "do_check", [](const reader_config_t &c) { return c.do_check(); },
          [](reader_config_t &c, bool v) { c.do_check(v); },
          "Do detector consistency check")
      .def_property(
          "verbose_check",
          [](const reader_config_t &c) { return c.verbose_check(); },
          [](reader_config_t &c, bool v) { c.verbose_check(v); },
          "Verbosity of the detector consistency check")
      .def(
          "add_file",
          [](reader_config_t &c, const std::string &f) -> reader_config_t & {
            return c.add_file(f);
          },
          py::arg("file_name"), py::return_value_policy::reference_internal,
          "Add an input file")
      .def("__repr__", &to_string<reader_config_t>);

  py::class_<intersection_config_t>(m, "IntersectionConfig")
      .def(py::init<>())
      .def_readwrite("min_mask_tolerance",
                     &intersection_config_t::min_mask_tolerance,
                     "Minimum mask tolerance")
      .def_readwrite("max_mask_tolerance",
                     &intersection_config_t::max_mask_tolerance,
                     "Maximum mask tolerance")
      .def_readwrite("mask_tolerance_scalor",
                     &intersection_config_t::mask_tolerance_scalor,
                     "Mask tolerance scale factor")
      .def_readwrite("path_tolerance", &intersection_config_t::path_tolerance,
                     "Tolerance to decide when a track is on a surface")
      .def_readwrite("overstep_tolerance",
                     &intersection_config_t::overstep_tolerance,
                     "How far behind the track position to look for candidates")
      .def("__repr__", &to_string<intersection_config_t>);

  py::class_<navigation_config_t>(m, "NavigationConfig")
      .def(py::init<>())
      .def_readwrite("intersection", &navigation_config_t::intersection,
                     "Intersection configuration")
      .def_readwrite("search_window", &navigation_config_t::search_window,
                     "Search window size for grid based acceleration structures")
      .def_readwrite("accumulated_error",
                     &navigation_config_t::accumulated_error,
                     "Percentage of total track path to assume as accumulated error")
      .def_readwrite("n_scattering_stddev",
                     &navigation_config_t::n_scattering_stddev,
                     "No. of standard deviations to assume to model the scattering noise")
      .def_readwrite("estimate_scattering_noise",
                     &navigation_config_t::estimate_scattering_noise,
                     "Add adaptive mask tolerance to navigation")
      .def("__repr__", &to_string<navigation_config_t>);

  py::class_<stepping_config_t>(m, "SteppingConfig")
      .def(py::init<>())
      .def_readwrite("min_stepsize", &stepping_config_t::min_stepsize,
                     "Minimum step size")
      .def_readwrite("rk_error_tol", &stepping_config_t::rk_error_tol,
                     "Runge-Kutta numeric error tolerance")
      .def_readwrite("step_constraint", &stepping_config_t::step_constraint,
                     "Step size constraint")
      .def_readwrite("path_limit", &stepping_config_t::path_limit,
                     "Maximum path length of track")
      .def_readwrite("max_rk_updates", &stepping_config_t::max_rk_updates,
                     "Maximum number of Runge-Kutta step trials")
      .def_readwrite("use_mean_loss", &stepping_config_t::use_mean_loss,
                     "Use mean energy loss (Bethe), otherwise use most probable energy loss (Landau)")
      .def_readwrite("use_eloss_gradient",
                     &stepping_config_t::use_eloss_gradient,
                     "Use energy loss gradient in error propagation")
      .def_readwrite("use_field_gradient",
                     &stepping_config_t::use_field_gradient,
                     "Use field gradient in error propagation")
      .def_readwrite("do_covariance_transport",
                     &stepping_config_t::do_covariance_transport,
                     "Do covariance transport")
      .def("__repr__", &to_string<stepping_config_t>);

  py::class_<propagation_config_t>(m, "PropagationConfig")
      .def(py::init<>())
      .def_readwrite("navigation", &propagation_config_t::navigation,
                     "Navigation configuration")
      .def_readwrite("stepping", &propagation_config_t::stepping,
                     "Stepping configuration")
      .def("__repr__", &to_string<propagation_config_t>);

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
      .def("__repr__",
           [](const detector_handle &d) { return to_string(d.detector); });

  m.def("readDetector", &read_detector, py::arg("config"),
        "Read a detector as configured by a DetectorReaderConfig");
}

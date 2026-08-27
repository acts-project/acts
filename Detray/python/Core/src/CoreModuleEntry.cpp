// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// Detray core include(s)
#include "detray/core/detector.hpp"

// Detray navigation include(s)
#include "detray/navigation/volume_graph.hpp"

// Detray propagation include(s)
#include "detray/propagator/propagation_config.hpp"

// Detray algebra plugin + detector metadata
#include "algebra/array.hpp"
#include "detray/definitions/algebra.hpp"
#include "detray/detectors/default_metadata.hpp"

// Vecmem include(s)
#include <vecmem/memory/host_memory_resource.hpp>
#include <vecmem/memory/memory_resource.hpp>

// Pybind11 include(s)
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

// System include(s)
#include <array>
#include <cstddef>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>

namespace py = pybind11;

#define STRINGIFY(x) #x
/// Expand macro and turn its value into a string.
#define STRINGIFY_HELPER(x) STRINGIFY(x)

namespace {

using intersection_config_t = detray::intersection::config;
using navigation_config_t = detray::navigation::config;
using stepping_config_t = detray::stepping::config;
using propagation_config_t = detray::propagation::config;
using scalar_t = DETRAY_CUSTOM_SCALARTYPE;

template <typename T>
std::string to_string(const T &obj) {
  std::ostringstream os;
  os << obj;
  return os.str();
}

/// Bind @p Container , a read-only vector with elements type @p Element.
template <typename Container, typename Element>
void bind_const_vector(py::module_ &m, const char *name) {
  py::class_<Container>(m, name)
      .def("__len__", [](const Container &c) { return c.size(); })
      .def(
          "__getitem__",
          [](const Container &c, py::ssize_t i) -> const Element & {
            const auto n = static_cast<py::ssize_t>(c.size());
            if (i < 0) {
              i += n;
            }
            if (i < 0 || i >= n) {
              throw py::index_error();
            }
            return c[static_cast<std::size_t>(i)];
          },
          py::arg("index"), py::return_value_policy::reference_internal,
          "Element at a given index")
      .def(
          "__iter__",
          [](const Container &c) {
            return py::make_iterator<
                py::return_value_policy::reference_internal,
                decltype(c.begin()), decltype(c.end()), const Element &>(
                c.begin(), c.end());
          },
          py::keep_alive<0, 1>(), "Iterate over the elements");
}

using algebra_t = detray::array<scalar_t>;
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
using volume_graph_t = detray::volume_graph<detector_t>;

}  // namespace

// Treat the detector's volume/portal containers as opaque so that returning
// them by reference exposes the live storage instead of copying it into a
// Python list (the default for std::vector via pybind11/stl.h).
PYBIND11_MAKE_OPAQUE(volume_container_t)
PYBIND11_MAKE_OPAQUE(surface_container_t)

PYBIND11_MODULE(DetrayPythonBindings, m) {
  m.doc() = "Detray core bindings";

  py::class_<intersection_config_t>(m, "IntersectionConfig")
      .def(py::init<>())
      .def_readwrite("minMaskTolerance",
                     &intersection_config_t::min_mask_tolerance,
                     "Minimum mask tolerance")
      .def_readwrite("maxMaskTolerance",
                     &intersection_config_t::max_mask_tolerance,
                     "Maximum mask tolerance")
      .def_readwrite("maskToleranceScalor",
                     &intersection_config_t::mask_tolerance_scalor,
                     "Mask tolerance scale factor")
      .def_readwrite("pathTolerance", &intersection_config_t::path_tolerance,
                     "Tolerance to decide when a track is on a surface")
      .def_readwrite("overstepTolerance",
                     &intersection_config_t::overstep_tolerance,
                     "How far behind the track position to look for candidates")
      .def("__repr__", &to_string<intersection_config_t>);

  py::class_<navigation_config_t>(m, "NavigationConfig")
      .def(py::init<>())
      .def_readwrite("intersection", &navigation_config_t::intersection,
                     "Intersection configuration")
      .def_readwrite(
          "searchWindow", &navigation_config_t::search_window,
          "Search window size for grid based acceleration structures")
      .def_readwrite(
          "accumulatedError", &navigation_config_t::accumulated_error,
          "Percentage of total track path to assume as accumulated error")
      .def_readwrite(
          "nScatteringStddev", &navigation_config_t::n_scattering_stddev,
          "No. of standard deviations to assume to model the scattering noise")
      .def_readwrite("estimateScatteringNoise",
                     &navigation_config_t::estimate_scattering_noise,
                     "Add adaptive mask tolerance to navigation")
      .def("__repr__", &to_string<navigation_config_t>);

  py::class_<stepping_config_t>(m, "SteppingConfig")
      .def(py::init<>())
      .def_readwrite("minStepsize", &stepping_config_t::min_stepsize,
                     "Minimum step size")
      .def_readwrite("rkErrorTol", &stepping_config_t::rk_error_tol,
                     "Runge-Kutta numeric error tolerance")
      .def_readwrite("stepConstraint", &stepping_config_t::step_constraint,
                     "Step size constraint")
      .def_readwrite("pathLimit", &stepping_config_t::path_limit,
                     "Maximum path length of track")
      .def_readwrite("maxRkUpdates", &stepping_config_t::max_rk_updates,
                     "Maximum number of Runge-Kutta step trials")
      .def_readwrite("useMeanLoss", &stepping_config_t::use_mean_loss,
                     "Use mean energy loss (Bethe), otherwise use most "
                     "probable energy loss (Landau)")
      .def_readwrite("useElossGradient", &stepping_config_t::use_eloss_gradient,
                     "Use energy loss gradient in error propagation")
      .def_readwrite("useFieldGradient", &stepping_config_t::use_field_gradient,
                     "Use field gradient in error propagation")
      .def_readwrite("doCovarianceTransport",
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
  bind_const_vector<surface_store_t, surface_descriptor_t>(m, "SurfaceStore");
  bind_const_vector<volume_container_t, volume_descriptor_t>(m,
                                                             "VolumeContainer");
  bind_const_vector<surface_container_t, surface_descriptor_t>(
      m, "SurfaceContainer");
  py::class_<geometry_context_t>(m, "GeometryContext")
      .def(py::init<>())
      .def(py::init<detray::dindex>(), py::arg("index"))
      .def_property_readonly("index", &geometry_context_t::get,
                             "Index into the geometry data store")
      .def("__repr__", [](const geometry_context_t &gctx) {
        return "GeometryContext(index=" + std::to_string(gctx.get()) + ")";
      });
  py::class_<transform_store_t>(m, "TransformStore");
  py::class_<mask_store_t>(m, "MaskStore");
  py::class_<material_store_t>(m, "MaterialStore");
  py::class_<accelerator_store_t>(m, "AcceleratorStore");
  py::class_<detray::name_map>(m, "NameMap")
      .def(py::init<>())
      .def_property(
          "detectorName",
          [](const detray::name_map &n) { return n.get_detector_name(); },
          [](detray::name_map &n, std::string_view name) {
            n.set_detector_name(name);
          },
          "Name of the detector")
      .def("empty", &detray::name_map::empty,
           "Whether no volume names are mapped")
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
          [](detray::name_map &n, detray::dindex index,
             const std::string &name) { n.emplace(index, name); },
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
      .def("clearNames", &detray::name_map::clear_names,
           "Clear volume names, keep the detector name")
      .def("__repr__", [](const detray::name_map &n) {
        return "NameMap(detectorName='" + n.get_detector_name() + "')";
      });

  py::class_<vecmem::memory_resource, std::shared_ptr<vecmem::memory_resource>>(
      m, "MemoryResource");
  py::class_<vecmem::host_memory_resource, vecmem::memory_resource,
             std::shared_ptr<vecmem::host_memory_resource>>(
      m, "HostMemoryResource")
      .def(py::init<>());

  py::class_<detector_t>(
      m, "DetectorDefaultMetadata" STRINGIFY_HELPER(DETRAY_CUSTOM_SCALARTYPE))
      .def(
          "name",
          [](const detector_t &d, const detray::name_map &names) {
            return d.name(names);
          },
          py::arg("names"), "Detector name")
      .def_property_readonly(
          "volumes",
          [](const detector_t &d) -> const volume_container_t & {
            return d.volumes();
          },
          py::return_value_policy::reference_internal, "All volumes")
      .def_property_readonly(
          "surfaces",
          [](const detector_t &d) -> const surface_store_t & {
            return d.surfaces();
          },
          py::return_value_policy::reference_internal, "All surfaces")
      .def_property_readonly(
          "portals",
          [](const detector_t &d) -> const surface_container_t & {
            return d.portals();
          },
          py::return_value_policy::reference_internal, "All portals")
      .def_property_readonly(
          "transformStore",
          [](const detector_t &d) -> const transform_store_t & {
            return d.transform_store();
          },
          py::return_value_policy::reference_internal, "Transform store")
      .def_property_readonly(
          "maskStore",
          [](const detector_t &d) -> const mask_store_t & {
            return d.mask_store();
          },
          py::return_value_policy::reference_internal, "Mask store")
      .def_property_readonly(
          "materialStore",
          [](const detector_t &d) -> const material_store_t & {
            return d.material_store();
          },
          py::return_value_policy::reference_internal, "Material store")
      .def_property_readonly(
          "acceleratorStore",
          [](const detector_t &d) -> const accelerator_store_t & {
            return d.accelerator_store();
          },
          py::return_value_policy::reference_internal, "Accelerator store")
      .def("__repr__", [](const detector_t &d) { return to_string(d); });

  py::class_<volume_graph_t>(m, "VolumeGraph")
      .def(py::init<const detector_t &>(), py::arg("detector"),
           py::keep_alive<1, 2>(),
           "Build the volume graph of a detector. The detector is kept alive "
           "for as long as the graph is used")
      .def("toDotString", &volume_graph_t::to_dot_string,
           "The volume linking description in DOT syntax")
      .def("__repr__", &volume_graph_t::to_string);
}

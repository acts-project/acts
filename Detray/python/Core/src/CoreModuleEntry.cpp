// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// Detray core include(s)
#include "detray/python/detector.hpp"
#include "detray/utils/consistency_checker.hpp"

// Detray navigation include(s)
#include "detray/navigation/volume_graph.hpp"

// Detray propagation include(s)
#include "detray/propagator/propagation_config.hpp"

// Detray algebra plugin + detector metadata
#include "algebra/array.hpp"
#include "detray/definitions/algebra.hpp"
#include "detray/detectors/default_metadata.hpp"
#include "detray/detectors/toy_metadata.hpp"

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

namespace {

using intersection_config_t = detray::intersection::config;
using navigation_config_t = detray::navigation::config;
using stepping_config_t = detray::stepping::config;
using propagation_config_t = detray::propagation::config;

template <typename T>
std::string to_string(const T& obj) {
  std::ostringstream os;
  os << obj;
  return os.str();
}

using geometry_context_t = detray::geometry_context;

struct SurfaceInfo {
  std::uint64_t identifier;
  detray::dindex index, volume, transformIndex;
  bool sensitive, portal;
};
struct VolumeInfo {
  detray::dindex index, transformIndex;
  std::array<detray::dindex, 2> surfaceRange;
};

/// Keep the Python detector alive, including its allocator lifetime attachment.
struct GeometryCollection {
  py::object owner;
  enum class Kind { volumes, surfaces, portals } kind;
  std::size_t size() const {
    return detray::python::visit(
        owner.cast<const detray::host_detector&>(),
        [this]<typename T>(const T& detector) -> std::size_t {
          if (kind == Kind::volumes) {
            return detector.volumes().size();
          }
          if (kind == Kind::portals) {
            return detector.portals().size();
          }
          return detector.surfaces().size();
        });
  }
  py::object at(py::ssize_t index) const {
    auto count = static_cast<py::ssize_t>(size());
    if (index < 0) {
      index += count;
    }
    if (index < 0 || index >= count) {
      throw py::index_error();
    }
    return detray::python::visit(
        owner.cast<const detray::host_detector&>(),
        [this, index]<typename T>(const T& detector) -> py::object {
          if (kind == Kind::volumes) {
            const auto& volume = detector.volumes()[index];
            const auto range = volume.full_sf_range();
            return py::cast(VolumeInfo{
                volume.index(), volume.transform(), {range[0], range[1]}});
          }
          auto describe = [](const auto& surface) {
            return py::cast(
                SurfaceInfo{surface.identifier().value(), surface.index(),
                            surface.volume(), surface.transform(),
                            surface.is_sensitive(), surface.is_portal()});
          };
          if (kind == Kind::portals) {
            return describe(detector.portals()[index]);
          }
          return describe(detector.surfaces()[index]);
        });
  }
};

struct VolumeGraph {
  std::string dot, description;
  explicit VolumeGraph(const detray::host_detector& detector) {
    detray::python::visit(detector, [this]<typename T>(const T& concrete) {
      detray::volume_graph<T> graph{concrete};
      dot = graph.to_dot_string();
      description = graph.to_string();
    });
  }
};

}  // namespace

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

  py::class_<SurfaceInfo>(m, "SurfaceInfo")
      .def_readonly("identifier", &SurfaceInfo::identifier)
      .def_readonly("index", &SurfaceInfo::index)
      .def_readonly("volume", &SurfaceInfo::volume)
      .def_readonly("transformIndex", &SurfaceInfo::transformIndex)
      .def_readonly("sensitive", &SurfaceInfo::sensitive)
      .def_readonly("portal", &SurfaceInfo::portal);
  py::class_<VolumeInfo>(m, "VolumeInfo")
      .def_readonly("index", &VolumeInfo::index)
      .def_readonly("transformIndex", &VolumeInfo::transformIndex)
      .def_readonly("surfaceRange", &VolumeInfo::surfaceRange);
  py::class_<GeometryCollection>(m, "GeometryCollection")
      .def("__len__", &GeometryCollection::size)
      .def("__getitem__", &GeometryCollection::at);

  py::class_<geometry_context_t>(m, "GeometryContext")
      .def(py::init<>())
      .def(py::init<detray::dindex>(), py::arg("index"))
      .def_property_readonly("index", &geometry_context_t::get,
                             "Index into the geometry data store")
      .def("__repr__", [](const geometry_context_t& gctx) {
        return "GeometryContext(index=" + std::to_string(gctx.get()) + ")";
      });
  py::class_<detray::name_map>(m, "NameMap")
      .def(py::init<>())
      .def_property(
          "detectorName",
          [](const detray::name_map& n) { return n.get_detector_name(); },
          [](detray::name_map& n, std::string_view name) {
            n.set_detector_name(name);
          },
          "Name of the detector")
      .def("empty", &detray::name_map::empty,
           "Whether no volume names are mapped")
      .def(
          "__contains__",
          [](const detray::name_map& n, detray::dindex index) {
            return n.contains(index);
          },
          py::arg("index"), "Whether a volume index is mapped")
      .def(
          "__contains__",
          [](const detray::name_map& n, std::string_view name) {
            return n.contains(name);
          },
          py::arg("name"), "Whether a volume name is mapped")
      .def(
          "__setitem__",
          [](detray::name_map& n, detray::dindex index,
             const std::string& name) { n.emplace(index, name); },
          py::arg("index"), py::arg("name"), "Map a volume index to a name")
      .def(
          "__getitem__",
          [](const detray::name_map& n, detray::dindex index) -> std::string {
            try {
              return n.at(index);
            } catch (const std::out_of_range&) {
              throw py::key_error(std::to_string(index));
            }
          },
          py::arg("index"), "Volume name at a volume index")
      .def(
          "__getitem__",
          [](const detray::name_map& n,
             std::string_view name) -> detray::dindex {
            try {
              return n.at(name);
            } catch (const std::out_of_range&) {
              throw py::key_error(std::string{name});
            }
          },
          py::arg("name"), "Volume index at a volume name")
      .def("clear", &detray::name_map::clear, "Clear detector and volume names")
      .def("clearNames", &detray::name_map::clear_names,
           "Clear volume names, keep the detector name")
      .def("__repr__", [](const detray::name_map& n) {
        return "NameMap(detectorName='" + n.get_detector_name() + "')";
      });

  py::class_<vecmem::memory_resource, std::shared_ptr<vecmem::memory_resource>>(
      m, "MemoryResource");
  py::class_<vecmem::host_memory_resource, vecmem::memory_resource,
             std::shared_ptr<vecmem::host_memory_resource>>(
      m, "HostMemoryResource")
      .def(py::init<>());

  py::class_<detray::host_detector, std::shared_ptr<detray::host_detector>>(
      m, "Detector")
      .def_property_readonly(
          "metadataId",
          [](const detray::host_detector& detector) {
            return detray::python::visit(detector, []<typename T>(const T&) {
              return std::string(detray::registry::metadata_id<
                                 detray::detector_metadata_t<T>>());
            });
          })
      .def(
          "name",
          [](const detray::host_detector& detector,
             const detray::name_map& names) {
            return detray::python::visit(detector,
                                         [&]<typename T>(const T& concrete) {
                                           return concrete.name(names);
                                         });
          },
          py::arg("names"))
      .def_property_readonly("volumes",
                             [](py::object self) {
                               return GeometryCollection{
                                   std::move(self),
                                   GeometryCollection::Kind::volumes};
                             })
      .def_property_readonly("surfaces",
                             [](py::object self) {
                               return GeometryCollection{
                                   std::move(self),
                                   GeometryCollection::Kind::surfaces};
                             })
      .def_property_readonly("portals",
                             [](py::object self) {
                               return GeometryCollection{
                                   std::move(self),
                                   GeometryCollection::Kind::portals};
                             })
      .def("checkConsistency",
           [](const detray::host_detector& detector) {
             detray::python::visit(detector, []<typename T>(const T& concrete) {
               detray::detail::check_consistency(concrete);
             });
           })
      .def("__repr__", [](const detray::host_detector& detector) {
        return detray::python::visit(
            detector,
            []<typename T>(const T& concrete) { return to_string(concrete); });
      });
  m.attr("availableMetadata") = py::cast(detray::registry::metadata_ids);
  py::class_<VolumeGraph>(m, "VolumeGraph")
      .def(py::init<const detray::host_detector&>(), py::arg("detector"))
      .def("toDotString", [](const VolumeGraph& graph) { return graph.dot; })
      .def("__repr__",
           [](const VolumeGraph& graph) { return graph.description; });
}

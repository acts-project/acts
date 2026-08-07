// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// Detray core include(s)
#include "detray/core/detector.hpp"

// Detray propagation include(s)
#include "detray/propagator/propagation_config.hpp"

// Detray event generator include(s)
#include "detray/test/common/event_generator/uniform_track_generator_config.hpp"

// Detray test include(s)
#include "detray/test/cpu/material_scan.hpp"

// Detray algebra plugin + detector metadata
#include "algebra/array.hpp"
#include "detray/definitions/algebra.hpp"
#include "detray/detectors/default_metadata.hpp"

// Pybind11 include(s)
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

// System include(s)
#include <cstddef>
#include <sstream>
#include <string>

namespace py = pybind11;

namespace {

using scalar_t = DETRAY_CUSTOM_SCALARTYPE;
using algebra_t = detray::array<scalar_t>;
using detector_t = detray::detector<detray::default_metadata<algebra_t>>;
using material_scan_config_t = detray::test::material_scan<detector_t>::config;
using propagation_config_t = detray::propagation::config;
using track_generator_config_t =
    detray::uniform_track_generator_config<scalar_t>;

template <typename T>
std::string to_string(const T &obj) {
  std::ostringstream os;
  os << obj;
  return os.str();
}

}  // namespace

PYBIND11_MODULE(DetrayTestsPythonBindings, m) {
  m.doc() = "Detray tests bindings";

  py::class_<material_scan_config_t>(m, "MaterialScanConfig")
      .def(py::init<>())
      .def_property(
          "name", [](const material_scan_config_t &c) { return c.name(); },
          [](material_scan_config_t &c, const std::string &n) { c.name(n); },
          "Name of the test")
      .def_property(
          "materialFile",
          [](const material_scan_config_t &c) { return c.material_file(); },
          [](material_scan_config_t &c, const std::string &f) {
            c.material_file(f);
          },
          "Name of the output file with the material traces")
      .def_property(
          "overlapsRemoval",
          [](const material_scan_config_t &c) { return c.overlaps_removal(); },
          [](material_scan_config_t &c, bool o) { c.overlaps_removal(o); },
          "Perform overlaps removal")
      .def_property(
          "overlapsTol",
          [](const material_scan_config_t &c) { return c.overlaps_tol(); },
          [](material_scan_config_t &c, scalar_t t) { c.overlaps_tol(t); },
          "Tolerance for considering surfaces to be overlapping")
      .def_property(
          "trackGenerator",
          [](material_scan_config_t &c) -> track_generator_config_t & {
            return c.track_generator();
          },
          [](material_scan_config_t &c, const track_generator_config_t &v) {
            c.track_generator() = v;
          },
          py::return_value_policy::reference_internal,
          "Track generator configuration")
      .def_property(
          "tol", [](const material_scan_config_t &c) { return c.tol(); },
          [](material_scan_config_t &c, scalar_t t) { c.tol(t); },
          "Tolerance to compare two floating point values")
      .def_property(
          "propagation",
          [](material_scan_config_t &c) -> propagation_config_t & {
            return c.propagation();
          },
          [](material_scan_config_t &c, const propagation_config_t &v) {
            c.propagation() = v;
          },
          py::return_value_policy::reference_internal,
          "Propagation configuration")
      .def("__repr__", &to_string<material_scan_config_t>);
}

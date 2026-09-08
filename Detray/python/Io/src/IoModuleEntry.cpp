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
#include <vecmem/memory/memory_resource.hpp>

// Pybind11 include(s)
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

// System include(s)
#include <memory>
#include <sstream>
#include <string>
#include <utility>

namespace py = pybind11;

namespace {

using reader_config_t = detray::io::detector_reader_config;
using scalar_t = DETRAY_CUSTOM_SCALARTYPE;
using algebra_t = detray::array<scalar_t>;
using detector_t = detray::detector<detray::default_metadata<algebra_t>>;

template <typename T>
std::string to_string(const T &obj) {
  std::ostringstream os;
  os << obj;
  return os.str();
}

/// Read a detector (default metadata) into @p mr as configured by @p cfg .
///
/// The memory resource object is automatically kept alive at least as long as
/// the returned detector object.
std::pair<py::object, detray::name_map> read_detector(
    std::shared_ptr<vecmem::memory_resource> mr, const reader_config_t &cfg) {
  auto [det, names] = detray::io::read_detector<detector_t>(*mr, cfg);

  py::object detector = py::cast(std::move(det));
  py::detail::keep_alive_impl(detector, py::cast(std::move(mr)));

  return {std::move(detector), std::move(names)};
}

}  // namespace

PYBIND11_MODULE(DetrayIoPythonBindings, m) {
  m.doc() = "Detray io bindings";

  py::class_<reader_config_t>(m, "DetectorReaderConfig")
      .def(py::init<>())
      .def_property_readonly("files", &reader_config_t::files, "Input files")
      .def_property(
          "doCheck", [](const reader_config_t &c) { return c.do_check(); },
          [](reader_config_t &c, bool v) { c.do_check(v); },
          "Do detector consistency check")
      .def_property(
          "verboseCheck",
          [](const reader_config_t &c) { return c.verbose_check(); },
          [](reader_config_t &c, bool v) { c.verbose_check(v); },
          "Verbosity of the detector consistency check")
      .def(
          "addFile",
          [](reader_config_t &c, const std::string &f) -> reader_config_t & {
            return c.add_file(f);
          },
          py::arg("fileName"), py::return_value_policy::reference_internal,
          "Add an input file")
      .def("__repr__", &to_string<reader_config_t>);

  m.def("readDetector", &read_detector, py::arg("memoryResource"),
        py::arg("config"),
        "Read a detector using the given memory resource as configured by a "
        "DetectorReaderConfig. The returned detector keeps the memory resource "
        "alive for as long as it is used");
}

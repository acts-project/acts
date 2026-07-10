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
#include <string>
#include <utility>

namespace py = pybind11;

namespace {

using algebra_t = detray::array<DETRAY_CUSTOM_SCALARTYPE>;
using detector_t = detray::detector<detray::default_metadata<algebra_t>>;
using name_map_t = detector_t::name_map;

/// Owns a detector together with the memory resource its data lives in.
struct detector_handle {
  std::unique_ptr<vecmem::host_memory_resource> memory_resource;
  detector_t detector;
};

/// Read a detector (default metadata) from a JSON file.
std::pair<detector_handle, name_map_t> read_detector(
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

  py::class_<detector_handle>(m, "Detector")
      .def(
          "n_volumes",
          [](const detector_handle &d) { return d.detector.volumes().size(); },
          "Number of volumes in the detector")
      .def(
          "n_surfaces",
          [](const detector_handle &d) { return d.detector.surfaces().size(); },
          "Number of surfaces in the detector");
  py::class_<name_map_t>(m, "NameMap");

  m.def("read_detector", &read_detector, py::arg("file_name"),
        "Read a detector from a JSON file");
}

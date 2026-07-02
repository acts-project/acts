// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// Detray IO include(s)
#include "detray/io/frontend/detector_writer.hpp"

// Detray test include(s)
#include "detray/test/common/build_toy_detector.hpp"
#include "detray/test/framework/types.hpp"

// Vecmem include(s)
#include <vecmem/memory/host_memory_resource.hpp>

// Pybind11 include(s)
#include <pybind11/pybind11.h>

// System include(s)
#include <stdexcept>
#include <string>

namespace py = pybind11;

namespace {

/// Build a toy detector and write it to JSON files in @param output_dir.
void generate_toy_detector(unsigned int barrel_layers,
                           unsigned int endcap_layers, bool material_use_maps,
                           const std::string& output_dir) {
  if (barrel_layers > 4) {
    throw std::invalid_argument(
        "Number of barrel layers out of range (0-4): got '" +
        std::to_string(barrel_layers) + "'");
  }
  if (endcap_layers > 7) {
    throw std::invalid_argument(
        "Number of endcap layers out of range (0-7): got '" +
        std::to_string(endcap_layers) + "'");
  }

  detray::toy_det_config<detray::test::scalar> cfg{};
  cfg.n_brl_layers(barrel_layers);
  cfg.n_edc_layers(endcap_layers);
  cfg.use_material_maps(material_use_maps);

  vecmem::host_memory_resource host_mr;
  auto [det, names] =
      detray::build_toy_detector<detray::test::algebra>(host_mr, cfg);

  detray::io::detector_writer_config writer_cfg{};
  writer_cfg.format(detray::io::format::json).replace_files(false);
  writer_cfg.path(output_dir);

  detray::io::write_detector(det, names, writer_cfg);
}

}  // namespace

PYBIND11_MODULE(DetrayExamplesPythonBindings, m) {
  m.doc() = "Detray example bindings";

  m.def("generate_toy_detector", &generate_toy_detector,
        py::arg("barrel_layers") = 4u, py::arg("endcap_layers") = 3u,
        py::arg("material_use_maps") = false,
        py::arg("output_dir") = "./toy_detector/",
        "Build a toy detector and write it to JSON files in output_dir");
}

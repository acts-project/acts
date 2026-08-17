// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// Detray plugin include(s)
#include "detray/plugins/svgtools/writer.hpp"

// Actsvg include(s)
#include <actsvg/core/draw.hpp>
#include <actsvg/core/style.hpp>
#include <actsvg/core/svg.hpp>

// Pybind11 include(s)
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

// System include(s)
#include <array>
#include <string>
#include <vector>

namespace py = pybind11;

PYBIND11_MODULE(DetraySvgtoolsPythonBindings, m) {
  m.doc() = "Detray svgtools bindings";

  // Svg primitives
  {
    // Duplicates actsvg's binding. Registered module-locally, so that importing
    // this module alongside actsvg does not clash over the same C++ type
    py::class_<actsvg::svg::object>(m, "SvgObject", py::module_local())
        .def_readwrite("id", &actsvg::svg::object::_id,
                       "Identification string");

    m.def(
        "writeSvg",
        [](const std::string &path,
           const std::vector<actsvg::svg::object> &svgs,
           bool replace) { detray::svgtools::write_svg(path, svgs, replace); },
        py::arg("path"), py::arg("svgs"), py::arg("replace") = true,
        "Write svg objects to '<path>.svg'. Their ids must be unique");

    m.def(
        "drawAxes",
        [](const std::string &name,
           const std::array<actsvg::scalar, 2> &x_range,
           const std::array<actsvg::scalar, 2> &y_range,
           const std::string &x_label, const std::string &y_label,
           unsigned int font_size) {
          actsvg::style::font font{};
          font._size = font_size;

          // The default stroke is black
          return actsvg::draw::x_y_axes(name, x_range, y_range,
                                        actsvg::style::stroke{}, x_label,
                                        y_label, font);
        },
        py::arg("name"), py::arg("xRange"), py::arg("yRange"),
        py::arg("xLabel") = "", py::arg("yLabel") = "",
        py::arg("fontSize") = 12u, "Draw a pair of x-y axes");
  }
}

// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Visualization/IVisualization3D.hpp"
#include "Acts/Visualization/ViewConfig.hpp"
#include "ActsPython/Utilities/Context.hpp"
#include "ActsPython/Utilities/Macros.hpp"
#include "ActsPython/Utilities/Patchers.hpp"

#include <pybind11/pybind11.h>

namespace py = pybind11;
using namespace Acts;

namespace ActsPython {
/// This adds the definitions from Core/Visualization to the python module
/// @param ctx the context container for the python modules
void addVisualization(Context& ctx) {
  auto& m = ctx.get("main");

  {
    auto c = py::class_<ViewConfig>(m, "ViewConfig").def(py::init<>());

    ACTS_PYTHON_STRUCT(c, visible, color, offset, lineThickness,
                       surfaceThickness, quarterSegments, triangulate,
                       outputName);

    patchKwargsConstructor(c);

    py::class_<Color>(m, "Color")
        .def(py::init<>())
        .def(py::init<int, int, int>())
        .def(py::init<double, double, double>())
        .def(py::init<std::string_view>())
        .def_readonly("rgb", &Color::rgb);
  }

  py::class_<IVisualization3D>(m, "IVisualization3D")
      .def("write", py::overload_cast<const std::filesystem::path&>(
                        &IVisualization3D::write, py::const_));
}
}  // namespace ActsPython
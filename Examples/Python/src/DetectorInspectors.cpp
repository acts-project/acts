// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/Plugins/Python/Utilities.hpp"
#include "ActsExamples/DetectorInspectors/SurfaceIndexing.hpp"

#include <memory>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;
using namespace ActsExamples;

namespace Acts::Python {
void addDetectorInspectors(Context& ctx) {
  auto [m, mex] = ctx.get("main", "examples");
  {
    auto ci = py::class_<CylindricalDetectorIndexing,
                         std::shared_ptr<CylindricalDetectorIndexing>>(
                  mex, "CylindricalDetectorIndexing")
                  .def(py::init<std::string>())
                  .def("inspect", &CylindricalDetectorIndexing::inspect);
  }
}
}  // namespace Acts::Python

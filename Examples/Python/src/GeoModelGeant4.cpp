// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Plugins/Python/Utilities.hpp"
#include "ActsExamples/GeoModelDetector/GeoModelGeant4DetectorConstruction.hpp"

#include <memory>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;
using namespace ActsExamples;

namespace Acts::Python {

void addGeoModelGeant4(Context& ctx) {
  auto [m, mex] = ctx.get("main", "examples");

  {
    py::class_<GeoModelGeant4DetectorConstruction,
               std::shared_ptr<GeoModelGeant4DetectorConstruction>>(
        m, "GeoModelGeant4DetectorConstruction");
  }
}

}  // namespace Acts::Python

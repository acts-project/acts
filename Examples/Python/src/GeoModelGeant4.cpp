// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Plugins/Python/Utilities.hpp"
#include "ActsExamples/GeoModelDetector/GeoModelGeant4DetectorConstructionFactory.hpp"

#include <memory>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;
using namespace ActsExamples;

namespace Acts::Python {

void addGeoModelGeant4(Context& ctx) {
  auto [m, mex] = ctx.get("main", "examples");

  {
    py::class_<GeoModelGeant4DetectorConstructionFactory,
               DetectorConstructionFactory,
               std::shared_ptr<GeoModelGeant4DetectorConstructionFactory>>(
        m, "GeoModelGeant4DetectorConstructionFactory")
        .def(py::init<const Acts::GeoModelTree&,
                      std::vector<std::shared_ptr<RegionCreator>>>(),
             py::arg("geoModelTree"),
             py::arg("regionCreators") =
                 std::vector<std::shared_ptr<RegionCreator>>());
  }
}

}  // namespace Acts::Python

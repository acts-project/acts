// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "ActsPlugins/Json/JsonMaterialDecorator.hpp"
#include "ActsPlugins/Json/JsonSurfacesReader.hpp"
#include "ActsPlugins/Json/MaterialMapJsonConverter.hpp"
#include "ActsPython/Utilities/Helpers.hpp"
#include "ActsPython/Utilities/Macros.hpp"

#include <memory>
#include <string>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;
using namespace pybind11::literals;

using namespace Acts;
using namespace ActsPython;
using namespace ActsExamples;

PYBIND11_MODULE(ActsPluginsPythonBindingsJson, json) {
  {
    py::class_<JsonMaterialDecorator, IMaterialDecorator,
               std::shared_ptr<JsonMaterialDecorator>>(json,
                                                       "JsonMaterialDecorator")
        .def(py::init<const MaterialMapJsonConverter::Config&,
                      const std::string&, Logging::Level, bool, bool>(),
             py::arg("rConfig"), py::arg("jFileName"), py::arg("level"),
             py::arg("clearSurfaceMaterial") = true,
             py::arg("clearVolumeMaterial") = true);
  }

  {
    auto cls =
        py::class_<MaterialMapJsonConverter>(json, "MaterialMapJsonConverter")
            .def(py::init<const MaterialMapJsonConverter::Config&,
                          Logging::Level>(),
                 py::arg("config"), py::arg("level"));

    auto c = py::class_<MaterialMapJsonConverter::Config>(cls, "Config")
                 .def(py::init<>());
    ACTS_PYTHON_STRUCT(c, context, processSensitives, processApproaches,
                       processRepresenting, processBoundaries, processVolumes,
                       processDenseVolumes, processNonMaterial);
  }

  {
    auto sjOptions =
        py::class_<JsonSurfacesReader::Options>(json, "SurfaceJsonOptions")
            .def(py::init<>());
    ACTS_PYTHON_STRUCT(sjOptions, inputFile, jsonEntryPath);

    json.def("readSurfaceHierarchyMapFromJson",
             JsonSurfacesReader::readHierarchyMap);

    json.def("readSurfaceVectorFromJson", JsonSurfacesReader::readVector);

    py::class_<JsonDetectorElement, DetectorElementBase,
               std::shared_ptr<JsonDetectorElement>>(json,
                                                     "JsonDetectorElement")
        .def("surface", [](JsonDetectorElement& self) {
          return self.surface().getSharedPtr();
        });

    json.def("readDetectorElementsFromJson",
             JsonSurfacesReader::readDetectorElements);
  }
}

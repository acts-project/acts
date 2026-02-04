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
#include "ActsPlugins/Json/SurfaceJsonConverter.hpp"
#include "ActsPython/Utilities/Helpers.hpp"
#include "ActsPython/Utilities/Macros.hpp"

#include <fstream>
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
    auto sjOptions = py::class_<SurfaceJsonConverter::Options>(
                         json, "SurfaceJsonConverterOptions")
                         .def(py::init<>());
    ACTS_PYTHON_STRUCT(sjOptions, transformOptions, writeMaterial,
                       writeVertices, portal);

    json.def("surfaceToJson", &SurfaceJsonConverter::toJson, py::arg("gctx"),
             py::arg("surface"),
             py::arg("options") = SurfaceJsonConverter::Options{});

    json.def(
        "writeSurfaceVector",
        [](const GeometryContext& gctx,
           const std::vector<std::shared_ptr<const Surface>>& surfaces,
           const SurfaceJsonConverter::Options& options,
           const std::string& fileName) {
          nlohmann::json jSurfaces = nlohmann::json::array();
          for (const auto& surface : surfaces) {
            jSurfaces.push_back(
                SurfaceJsonConverter::toJson(gctx, *surface, options));
          }
          std::ofstream out(fileName);
          out << jSurfaces.dump(2);
        },
        py::arg("gctx"), py::arg("surfaces"), py::arg("options"),
        py::arg("fileName"));
  }

  {
    auto srjOptions = py::class_<JsonSurfacesReader::Options>(
                          json, "SurfaceReaderJsonOptions")
                          .def(py::init<>());
    ACTS_PYTHON_STRUCT(srjOptions, inputFile, jsonEntryPath);

    json.def("readSurfaceHierarchyMap", JsonSurfacesReader::readHierarchyMap);

    json.def("readSurfaceVector", JsonSurfacesReader::readVector);

    py::class_<JsonDetectorElement, SurfacePlacementBase,
               std::shared_ptr<JsonDetectorElement>>(json,
                                                     "JsonDetectorElement")
        .def("surface", [](JsonDetectorElement& self) {
          return self.surface().getSharedPtr();
        });

    json.def("readDetectorElements", JsonSurfacesReader::readDetectorElements);
  }
}

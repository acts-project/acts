// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "ActsPlugins/Json/JsonMaterialDecorator.hpp"
#include "ActsPlugins/Json/JsonSurfacesReader.hpp"
#include "ActsPlugins/Json/MaterialMapJsonConverter.hpp"
#include "ActsPlugins/Json/TrackingGeometryJsonConverter.hpp"
#include "ActsPython/Utilities/Helpers.hpp"
#include "ActsPython/Utilities/Macros.hpp"

#include <memory>
#include <string>

#include <nlohmann/json.hpp>
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
                      const std::string&, Logging::Level>(),
             py::arg("rConfig"), py::arg("jFileName"), py::arg("level"));
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

    py::class_<JsonDetectorElement, SurfacePlacementBase,
               std::shared_ptr<JsonDetectorElement>>(json,
                                                     "JsonDetectorElement")
        .def("surface", [](JsonDetectorElement& self) {
          return self.surface().getSharedPtr();
        });

    json.def("readDetectorElementsFromJson",
             JsonSurfacesReader::readDetectorElements);
  }

  {
    auto cls = py::class_<TrackingGeometryJsonConverter>(
        json, "TrackingGeometryJsonConverter");

    py::class_<TrackingGeometryJsonConverter::Config>(cls, "Config")
        .def(py::init(&TrackingGeometryJsonConverter::Config::defaultConfig))
        .def_static("defaultConfig",
                    &TrackingGeometryJsonConverter::Config::defaultConfig);

    cls.def(py::init([](TrackingGeometryJsonConverter::Config config,
                        Acts::Logging::Level level) {
              return std::make_unique<TrackingGeometryJsonConverter>(
                  std::move(config),
                  Acts::getDefaultLogger("TrackingGeometryJsonConverter",
                                         level));
            }),
            py::arg("config") =
                TrackingGeometryJsonConverter::Config::defaultConfig(),
            py::arg("level") = Acts::Logging::INFO)
        .def(py::init([](TrackingGeometryJsonConverter::Config config,
                         std::unique_ptr<const Acts::Logger> logger) {
               return std::make_unique<TrackingGeometryJsonConverter>(
                   std::move(config), std::move(logger));
             }),
             py::arg("config") =
                 TrackingGeometryJsonConverter::Config::defaultConfig(),
             py::arg("logger"))
        .def(
            "toJson",
            [](const TrackingGeometryJsonConverter& self,
               const GeometryContext& gctx, const TrackingGeometry& geometry) {
              return self.toJson(gctx, geometry).dump();
            },
            py::arg("gctx"), py::arg("geometry"))
        .def(
            "fromJson",
            [](const TrackingGeometryJsonConverter& self,
               const GeometryContext& gctx, const std::string& encoded) {
              return self.fromJson(gctx, nlohmann::json::parse(encoded));
            },
            py::arg("gctx"), py::arg("encoded"));
  }
}

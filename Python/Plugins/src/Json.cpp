// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Detector/Detector.hpp"
#include "Acts/Detector/ProtoDetector.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "ActsPlugins/Json/DetectorJsonConverter.hpp"
#include "ActsPlugins/Json/JsonMaterialDecorator.hpp"
#include "ActsPlugins/Json/JsonSurfacesReader.hpp"
#include "ActsPlugins/Json/MaterialMapJsonConverter.hpp"
#include "ActsPlugins/Json/ProtoDetectorJsonConverter.hpp"
#include "ActsPython/Utilities/Helpers.hpp"
#include "ActsPython/Utilities/Macros.hpp"

#include <memory>
#include <string>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;
using namespace pybind11::literals;

PYBIND11_MODULE(ActsPluginsPythonBindingsJson, json) {
  using namespace Acts;
  using namespace ActsPython;

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

  // Gen2 Geometry code - to be removed when Gen3 is default
  using namespace Acts::Experimental;

  {
    py::class_<ProtoDetector>(json, "ProtoDetector")
        .def(py::init<>([](std::string pathName) {
          nlohmann::json jDetector;
          auto in = std::ifstream(pathName, std::ifstream::in);
          if (in.good()) {
            in >> jDetector;
            in.close();
          }
          ProtoDetector pDetector = jDetector["detector"];
          return pDetector;
        }));
  }

  {
    json.def("writeDetectorToJson",
             [](const GeometryContext& gctx, const Detector& detector,
                const std::string& name) -> void {
               auto jDetector = DetectorJsonConverter::toJson(gctx, detector);
               std::ofstream out;
               out.open(name + ".json");
               out << jDetector.dump(4);
               out.close();
             });
  }

  {
    json.def("writeDetectorToJsonDetray",
             [](const GeometryContext& gctx, const Detector& detector,
                const std::string& name) -> void {
               // Detray format test - manipulate for detray
               DetectorVolumeJsonConverter::Options detrayOptions;
               detrayOptions.transformOptions.writeIdentity = true;
               detrayOptions.transformOptions.transpose = true;
               detrayOptions.surfaceOptions.transformOptions =
                   detrayOptions.transformOptions;
               detrayOptions.portalOptions.surfaceOptions =
                   detrayOptions.surfaceOptions;

               auto jDetector = DetectorJsonConverter::toJsonDetray(
                   gctx, detector,
                   DetectorJsonConverter::Options{detrayOptions});

               // Write out the geometry, surface_grid, material
               auto jGeometry = jDetector["geometry"];
               auto jSurfaceGrids = jDetector["surface_grids"];
               auto jMaterial = jDetector["material"];

               std::ofstream out;
               out.open(name + "_geometry_detray.json");
               out << jGeometry.dump(4);
               out.close();

               out.open(name + "_surface_grids_detray.json");
               out << jSurfaceGrids.dump(4);
               out.close();

               out.open(name + "_material_detray.json");
               out << jMaterial.dump(4);
               out.close();
             });
  }

  {
    json.def(
        "readDetectorFromJson",
        [](const GeometryContext& gctx, const std::string& fileName) -> auto {
          auto in = std::ifstream(fileName,
                                  std::ifstream::in | std::ifstream::binary);
          nlohmann::json jDetectorIn;
          in >> jDetectorIn;
          in.close();

          return DetectorJsonConverter::fromJson(gctx, jDetectorIn);
        });
  }
}

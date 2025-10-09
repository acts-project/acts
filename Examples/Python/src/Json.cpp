// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Detector/Detector.hpp"
#include "Acts/Detector/ProtoDetector.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/Framework/ProcessCode.hpp"
#include "ActsExamples/Io/Json/JsonMaterialWriter.hpp"
#include "ActsExamples/Io/Json/JsonSurfacesWriter.hpp"
#include "ActsExamples/Io/Json/JsonTrackParamsLookupReader.hpp"
#include "ActsExamples/Io/Json/JsonTrackParamsLookupWriter.hpp"
#include "ActsPlugins/Json/DetectorJsonConverter.hpp"
#include "ActsPlugins/Json/JsonMaterialDecorator.hpp"
#include "ActsPlugins/Json/JsonSurfacesReader.hpp"
#include "ActsPlugins/Json/MaterialMapJsonConverter.hpp"
#include "ActsPlugins/Json/ProtoDetectorJsonConverter.hpp"
#include "ActsPython/Utilities/Helpers.hpp"
#include "ActsPython/Utilities/Macros.hpp"

#include <fstream>
#include <initializer_list>
#include <memory>
#include <string>
#include <tuple>
#include <vector>

#include <nlohmann/json.hpp>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace Acts {
class IMaterialDecorator;
}  // namespace Acts
namespace ActsExamples {
class IMaterialWriter;
class IWriter;

namespace Experimental {
class ITrackParamsLookupWriter;
}  // namespace Experimental

}  // namespace ActsExamples

namespace py = pybind11;
using namespace pybind11::literals;

using namespace Acts;
using namespace Acts::Experimental;
using namespace ActsExamples;

namespace ActsPython {
void addJson(Context& ctx) {
  auto [m, mex] = ctx.get("main", "examples");

  {
    py::class_<JsonMaterialDecorator, IMaterialDecorator,
               std::shared_ptr<JsonMaterialDecorator>>(m,
                                                       "JsonMaterialDecorator")
        .def(py::init<const MaterialMapJsonConverter::Config&,
                      const std::string&, Logging::Level, bool, bool>(),
             py::arg("rConfig"), py::arg("jFileName"), py::arg("level"),
             py::arg("clearSurfaceMaterial") = true,
             py::arg("clearVolumeMaterial") = true);
  }

  {
    auto cls =
        py::class_<MaterialMapJsonConverter>(m, "MaterialMapJsonConverter")
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
    py::enum_<JsonFormat>(mex, "JsonFormat")
        .value("NoOutput", JsonFormat::NoOutput)
        .value("Json", JsonFormat::Json)
        .value("Cbor", JsonFormat::Cbor)
        .value("All", JsonFormat::All);
  }

  {
    auto cls =
        py::class_<JsonMaterialWriter, IMaterialWriter,
                   std::shared_ptr<JsonMaterialWriter>>(mex,
                                                        "JsonMaterialWriter")
            .def(py::init<const JsonMaterialWriter::Config&, Logging::Level>(),
                 py::arg("config"), py::arg("level"))
            .def("writeMaterial", &JsonMaterialWriter::writeMaterial)
            .def("write", &JsonMaterialWriter::write)
            .def_property_readonly("config", &JsonMaterialWriter::config);

    auto c =
        py::class_<JsonMaterialWriter::Config>(cls, "Config").def(py::init<>());
    ACTS_PYTHON_STRUCT(c, converterCfg, fileName, writeFormat);
  }

  {
    using IWriter = ITrackParamsLookupWriter;
    using Writer = JsonTrackParamsLookupWriter;
    using Config = Writer::Config;

    auto cls = py::class_<Writer, IWriter, std::shared_ptr<Writer>>(
                   mex, "JsonTrackParamsLookupWriter")
                   .def(py::init<const Config&>(), py::arg("config"))
                   .def("writeLookup", &Writer::writeLookup)
                   .def_property_readonly("config", &Writer::config);

    auto c = py::class_<Config>(cls, "Config")
                 .def(py::init<>())
                 .def(py::init<const std::string&>(), py::arg("path"));
    ACTS_PYTHON_STRUCT(c, path);
  }

  {
    using IReader = ITrackParamsLookupReader;
    using Reader = JsonTrackParamsLookupReader;
    using Config = Reader::Config;

    auto cls = py::class_<Reader, IReader, std::shared_ptr<Reader>>(
                   mex, "JsonTrackParamsLookupReader")
                   .def(py::init<const Config&>(), py::arg("config"))
                   .def("readLookup", &Reader::readLookup)
                   .def_property_readonly("config", &Reader::config);

    auto c =
        py::class_<Config>(cls, "Config")
            .def(py::init<>())
            .def(
                py::init<std::unordered_map<GeometryIdentifier, const Surface*>,
                         std::pair<double, double>>(),
                py::arg("refLayers"), py::arg("bins"));
    ACTS_PYTHON_STRUCT(c, refLayers, bins);
  }

  {
    auto cls =
        py::class_<JsonSurfacesWriter, IWriter,
                   std::shared_ptr<JsonSurfacesWriter>>(mex,
                                                        "JsonSurfacesWriter")
            .def(py::init<const JsonSurfacesWriter::Config&, Logging::Level>(),
                 py::arg("config"), py::arg("level"))
            .def("write", &JsonSurfacesWriter::write)
            .def_property_readonly("config", &JsonSurfacesWriter::config);

    auto c =
        py::class_<JsonSurfacesWriter::Config>(cls, "Config").def(py::init<>());

    ACTS_PYTHON_STRUCT(c, trackingGeometry, outputDir, outputPrecision,
                       writeLayer, writeApproach, writeSensitive, writeBoundary,
                       writePerEvent, writeOnlyNames);
  }

  {
    py::class_<ProtoDetector>(mex, "ProtoDetector")
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
    auto sjOptions =
        py::class_<JsonSurfacesReader::Options>(m, "SurfaceJsonOptions")
            .def(py::init<>());
    ACTS_PYTHON_STRUCT(sjOptions, inputFile, jsonEntryPath);

    m.def("readSurfaceHierarchyMapFromJson",
          JsonSurfacesReader::readHierarchyMap);

    m.def("readSurfaceVectorFromJson", JsonSurfacesReader::readVector);

    py::class_<JsonDetectorElement, DetectorElementBase,
               std::shared_ptr<JsonDetectorElement>>(m, "JsonDetectorElement")
        .def("surface", [](JsonDetectorElement& self) {
          return self.surface().getSharedPtr();
        });

    m.def("readDetectorElementsFromJson",
          JsonSurfacesReader::readDetectorElements);
  }

  {
    mex.def("writeDetectorToJson",
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
    mex.def("writeDetectorToJsonDetray",
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
    mex.def(
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
}  // namespace ActsPython

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
#include "Acts/Plugins/Json/DetectorJsonConverter.hpp"
#include "Acts/Plugins/Json/JsonMaterialDecorator.hpp"
#include "Acts/Plugins/Json/JsonSurfacesReader.hpp"
#include "Acts/Plugins/Json/MaterialMapJsonConverter.hpp"
#include "Acts/Plugins/Json/ProtoDetectorJsonConverter.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/Framework/ProcessCode.hpp"
#include "ActsExamples/Io/Json/JsonDigitizationConfig.hpp"
#include "ActsExamples/Io/Json/JsonMaterialWriter.hpp"
#include "ActsExamples/Io/Json/JsonSurfacesWriter.hpp"
#include "ActsExamples/Io/Json/JsonTrackParamsLookupReader.hpp"
#include "ActsExamples/Io/Json/JsonTrackParamsLookupWriter.hpp"
#include "ActsExamples/Io/Json/JsonGeometryList.hpp"
#include "ActsPython/Utilities/Context.hpp"
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
using namespace ActsExamples;

namespace ActsPython {
void addJsonInputOutput(Context& ctx) {
  auto& mex = ctx.get("examples");

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
            .def(py::init<const JsonMaterialWriter::Config&,
                          Acts::Logging::Level>(),
                 py::arg("config"), py::arg("level"))
            .def("writeMaterial", &JsonMaterialWriter::writeMaterial)
            .def("write", &JsonMaterialWriter::write)
            .def_property_readonly("config", &JsonMaterialWriter::config);

    auto c =
        py::class_<JsonMaterialWriter::Config>(cls, "Config").def(py::init<>());
    ACTS_PYTHON_STRUCT(c, converterCfg, fileName, writeFormat);
  }

  {
    using IWriter = ActsExamples::ITrackParamsLookupWriter;
    using Writer = ActsExamples::JsonTrackParamsLookupWriter;
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
    using IReader = ActsExamples::ITrackParamsLookupReader;
    using Reader = ActsExamples::JsonTrackParamsLookupReader;
    using Config = Reader::Config;

    auto cls = py::class_<Reader, IReader, std::shared_ptr<Reader>>(
                   mex, "JsonTrackParamsLookupReader")
                   .def(py::init<const Config&>(), py::arg("config"))
                   .def("readLookup", &Reader::readLookup)
                   .def_property_readonly("config", &Reader::config);

    auto c = py::class_<Config>(cls, "Config")
                 .def(py::init<>())
                 .def(py::init<std::unordered_map<Acts::GeometryIdentifier,
                                                  const Acts::Surface*>,
                               std::pair<double, double>>(),
                      py::arg("refLayers"), py::arg("bins"));
    ACTS_PYTHON_STRUCT(c, refLayers, bins);
  }

  {
    auto cls =
        py::class_<JsonSurfacesWriter, IWriter,
                   std::shared_ptr<JsonSurfacesWriter>>(mex,
                                                        "JsonSurfacesWriter")
            .def(py::init<const JsonSurfacesWriter::Config&,
                          Acts::Logging::Level>(),
                 py::arg("config"), py::arg("level"))
            .def("write", &JsonSurfacesWriter::write)
            .def_property_readonly("config", &JsonSurfacesWriter::config);

    auto c =
        py::class_<JsonSurfacesWriter::Config>(cls, "Config").def(py::init<>());

    ACTS_PYTHON_STRUCT(c, trackingGeometry, outputDir, outputPrecision,
                       writeLayer, writeApproach, writeSensitive, writeBoundary,
                       writePerEvent, writeOnlyNames);
  }

  mex.def("readDigiConfigFromJson", ActsExamples::readDigiConfigFromJson);
  mex.def("writeDigiConfigToJson", ActsExamples::writeDigiConfigToJson);

  mex.def("readJsonGeometryList", ActsExamples::readJsonGeometryList);

}
}  // namespace ActsPython

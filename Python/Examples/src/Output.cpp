// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/Io/Csv/CsvBFieldWriter.hpp"
#include "ActsExamples/Io/Csv/CsvGnnGraphWriter.hpp"
#include "ActsExamples/Io/Csv/CsvMeasurementWriter.hpp"
#include "ActsExamples/Io/Csv/CsvParticleWriter.hpp"
#include "ActsExamples/Io/Csv/CsvProtoTrackWriter.hpp"
#include "ActsExamples/Io/Csv/CsvSeedWriter.hpp"
#include "ActsExamples/Io/Csv/CsvSimHitWriter.hpp"
#include "ActsExamples/Io/Csv/CsvSpacePointWriter.hpp"
#include "ActsExamples/Io/Csv/CsvSpacePointsBucketWriter.hpp"
#include "ActsExamples/Io/Csv/CsvTrackParameterWriter.hpp"
#include "ActsExamples/Io/Csv/CsvTrackWriter.hpp"
#include "ActsExamples/Io/Csv/CsvTrackingGeometryWriter.hpp"
#include "ActsExamples/Io/Csv/CsvVertexWriter.hpp"
#include "ActsExamples/Io/Obj/ObjPropagationStepsWriter.hpp"
#include "ActsExamples/Io/Obj/ObjSimHitWriter.hpp"
#include "ActsExamples/Io/Obj/ObjTrackingGeometryWriter.hpp"
#include "ActsExamples/MaterialMapping/IMaterialWriter.hpp"
#include "ActsExamples/TrackFinding/ITrackParamsLookupWriter.hpp"
#include "ActsPython/Utilities/Macros.hpp"

#include <memory>
#include <string>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/stl/filesystem.h>

namespace py = pybind11;
using namespace pybind11::literals;

using namespace Acts;
using namespace ActsExamples;

namespace {
template <CsvBFieldWriter::CoordinateType CType, bool Grid>
void register_csv_bfield_writer_binding(pybind11::class_<CsvBFieldWriter>& w) {
  std::string name =
      std::string(CType == CsvBFieldWriter::CoordinateType::XYZ ? "Xyz"
                                                                : "Rz") +
      std::string(Grid ? "Grid" : "Gridless");

  using Config = CsvBFieldWriter::Config<CType, Grid>;
  w.def_static((std::string("run") + name).c_str(),
               [](const Config& config, Logging::Level level) {
                 CsvBFieldWriter::run(config, level);
               },
               py::arg("config"), py::arg("level"));
  auto c = py::class_<Config>(w, (std::string("Config") + name).c_str())
               .def(py::init<>());
  ACTS_PYTHON_STRUCT(c, fileName, bField, range, bins);
}
}  // namespace

namespace ActsPython {

void addOutput(py::module& mex) {
  {
    using Writer = ObjTrackingGeometryWriter;
    auto w =
        py::class_<Writer, std::shared_ptr<Writer>>(mex,
                                                    "ObjTrackingGeometryWriter")
            .def(py::init<const Writer::Config&, Logging::Level>(),
                 py::arg("config"), py::arg("level"))
            .def("write",
                 py::overload_cast<const AlgorithmContext&,
                                   const TrackingGeometry&>(&Writer::write));

    auto c = py::class_<Writer::Config>(w, "Config").def(py::init<>());

    ACTS_PYTHON_STRUCT(c, outputScalor, outputPrecision, outputDir,
                       containerView, volumeView, sensitiveView, passiveView,
                       gridView);
  }

  ACTS_PYTHON_DECLARE_WRITER(ObjPropagationStepsWriter, mex,
                             "ObjPropagationStepsWriter", collection, outputDir,
                             outputScalor, outputPrecision);

  ACTS_PYTHON_DECLARE_WRITER(
      ObjSimHitWriter, mex, "ObjSimHitWriter", inputSimHits, outputDir,
      outputStem, outputPrecision, drawConnections, momentumThreshold,
      momentumThresholdTraj, nInterpolatedPoints, keepOriginalHits);

  ACTS_PYTHON_DECLARE_WRITER(CsvParticleWriter, mex, "CsvParticleWriter",
                             inputParticles, outputDir, outputStem,
                             outputPrecision);

  ACTS_PYTHON_DECLARE_WRITER(
      CsvMeasurementWriter, mex, "CsvMeasurementWriter", inputMeasurements,
      inputClusters, inputMeasurementSimHitsMap, outputDir, outputPrecision);

  ACTS_PYTHON_DECLARE_WRITER(CsvSimHitWriter, mex, "CsvSimHitWriter",
                             inputSimHits, outputDir, outputStem,
                             outputPrecision);

  ACTS_PYTHON_DECLARE_WRITER(CsvSpacePointWriter, mex, "CsvSpacePointWriter",
                             inputSpacepoints, outputDir, outputPrecision);

  ACTS_PYTHON_DECLARE_WRITER(CsvSpacePointsBucketWriter, mex,
                             "CsvSpacePointsBucketWriter", inputBuckets,
                             outputDir, outputPrecision);

  ACTS_PYTHON_DECLARE_WRITER(CsvTrackWriter, mex, "CsvTrackWriter", inputTracks,
                             outputDir, fileName, inputMeasurementParticlesMap,
                             outputPrecision, nMeasurementsMin,
                             truthMatchProbMin, ptMin);

  ACTS_PYTHON_DECLARE_WRITER(CsvSeedWriter, mex, "CsvSeedWriter",
                             inputTrackParameters, inputSimSeeds, inputSimHits,
                             inputMeasurementParticlesMap,
                             inputMeasurementSimHitsMap, fileName, outputDir);

  ACTS_PYTHON_DECLARE_WRITER(CsvVertexWriter, mex, "CsvVertexWriter",
                             inputVertices, outputDir, outputStem,
                             outputPrecision);

  ACTS_PYTHON_DECLARE_WRITER(
      CsvTrackingGeometryWriter, mex, "CsvTrackingGeometryWriter",
      trackingGeometry, outputDir, outputPrecision, writeSensitive,
      writeBoundary, writeSurfaceGrid, writeLayerVolume, writePerEvent);

  ACTS_PYTHON_DECLARE_WRITER(CsvTrackParameterWriter, mex,
                             "CsvTrackParameterWriter", inputTracks, outputDir,
                             outputStem, outputPrecision);

  ACTS_PYTHON_DECLARE_WRITER(CsvProtoTrackWriter, mex, "CsvProtoTrackWriter",
                             inputSpacepoints, inputPrototracks, outputDir);

  {
    using Writer = CsvBFieldWriter;

    auto w = py::class_<Writer>(mex, "CsvBFieldWriter");

    py::enum_<Writer::CoordinateType>(w, "CoordinateType")
        .value("rz", Writer::CoordinateType::RZ)
        .value("xyz", Writer::CoordinateType::XYZ);

    register_csv_bfield_writer_binding<Writer::CoordinateType::XYZ, true>(w);
    register_csv_bfield_writer_binding<Writer::CoordinateType::XYZ, false>(w);
    register_csv_bfield_writer_binding<Writer::CoordinateType::RZ, true>(w);
    register_csv_bfield_writer_binding<Writer::CoordinateType::RZ, false>(w);
  }

  ACTS_PYTHON_DECLARE_WRITER(CsvGnnGraphWriter, mex, "CsvGnnGraphWriter",
                             inputGraph, outputDir, outputStem);

  py::class_<IMaterialWriter, std::shared_ptr<IMaterialWriter>>(
      mex, "IMaterialWriter");

  py::class_<ITrackParamsLookupWriter,
             std::shared_ptr<ITrackParamsLookupWriter>>(
      mex, "ITrackParamsLookupWriter");
}

}  // namespace ActsPython

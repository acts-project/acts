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
#include "ActsExamples/Traccc/TracccBackend.hpp"
#include "ActsExamples/Traccc/TracccCopyAlgs.hpp"
#include "ActsExamples/Traccc/TracccCopyTool.hpp"
#include "ActsExamples/Traccc/TracccDetector.hpp"
#include "ActsExamples/Traccc/TracccEventData.hpp"
#include "ActsExamples/Traccc/TracccMagneticField.hpp"
#include "ActsPlugins/Covfie/FieldConversion.hpp"
#include "ActsPlugins/Detray/DetrayGeometryConverter.hpp"
#include "ActsPlugins/Detray/DetraySurfaceIdMap.hpp"
#include "ActsPython/Utilities/Helpers.hpp"
#include "ActsPython/Utilities/Macros.hpp"

#ifdef ACTS_EXAMPLES_TRACCC_CUDA
#include "ActsExamples/Traccc/Cuda/CudaTracccBackend.hpp"
#endif

#include <filesystem>
#include <memory>
#include <string>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/stl/filesystem.h>
#include <traccc/bfield/magnetic_field.hpp>
#include <traccc/geometry/detector.hpp>
#include <traccc/geometry/host_detector.hpp>
#include <traccc/io/read_magnetic_field.hpp>
#include <vecmem/memory/memory_resource.hpp>

namespace py = pybind11;
using namespace pybind11::literals;

using namespace ActsExamples;
using namespace ActsPython;

namespace {

void addBackends(py::module_& m) {
  py::class_<TracccBackend, std::shared_ptr<TracccBackend>>(m, "TracccBackend")
      .def_property_readonly("name", &TracccBackend::name);

#ifdef ACTS_EXAMPLES_TRACCC_CUDA
  py::class_<CudaTracccBackend, TracccBackend,
             std::shared_ptr<CudaTracccBackend>>(m, "CudaTracccBackend")
      .def(py::init([](int device, Acts::Logging::Level level) {
             CudaTracccBackend::Config cfg;
             cfg.device = device;
             return std::make_shared<CudaTracccBackend>(
                 cfg, Acts::getDefaultLogger("CudaTracccBackend", level));
           }),
           "device"_a = 0, "level"_a = Acts::Logging::INFO);
#endif
}

void addDetector(py::module_& m) {
  py::class_<traccc::host_detector, std::shared_ptr<traccc::host_detector>>(
      m, "TracccHostDetector");

  // Convert an ACTS tracking geometry into a type erased detray detector,
  // returned together with the ACTS <-> detray surface identifier map of the
  // conversion. The detray detector is move only, so it is moved into the
  // type erased container right after the conversion.
  m.def(
      "convertToTracccHostDetector",
      [](const ActsPlugins::DetrayGeometryConverter& converter,
         vecmem::memory_resource& mr, const Acts::GeometryContext& gctx,
         const std::shared_ptr<const Acts::TrackingGeometry>& trackingGeometry,
         const std::string& detectorName) {
        auto geometry = converter.convert<traccc::odd_detector::metadata>(
            mr, gctx, trackingGeometry, detectorName);
        auto hostDetector = std::make_shared<traccc::host_detector>();
        hostDetector->set<traccc::odd_detector>(std::move(*geometry.detector));
        auto surfaceIds = std::make_shared<ActsPlugins::DetraySurfaceIdMap>(
            std::move(geometry.surfaceIds));
        return std::make_pair(hostDetector, surfaceIds);
      },
      "converter"_a, "mr"_a, "gctx"_a, "trackingGeometry"_a,
      "detectorName"_a = "");

  auto d = py::class_<TracccDetector, std::shared_ptr<TracccDetector>>(
               m, "TracccDetector")
               .def(py::init([](const TracccDetector::Config& cfg,
                                Acts::Logging::Level level) {
                      return std::make_shared<TracccDetector>(
                          cfg, Acts::getDefaultLogger("TracccDetector", level));
                    }),
                    "config"_a, "level"_a = Acts::Logging::INFO)
               .def_property_readonly("config", &TracccDetector::config);

  auto c = py::class_<TracccDetector::Config>(d, "Config").def(py::init<>());
  ACTS_PYTHON_STRUCT(c, backend, trackingGeometry, detrayMemoryResource,
                     detrayDetector, surfaceIds);
}

void addMagneticField(py::module_& m) {
  py::class_<traccc::magnetic_field, std::shared_ptr<traccc::magnetic_field>>(
      m, "TracccHostMagneticField");

  // Build the traccc host field from an existing covfie field, or read it
  // from a traccc (covfie) field file
  m.def(
      "makeTracccHostMagneticField",
      [](const ActsPlugins::Covfie::ConstantField& field) {
        return std::make_shared<traccc::magnetic_field>(
            ActsPlugins::Covfie::ConstantField(field));
      },
      "field"_a);
  m.def(
      "readTracccHostMagneticField",
      [](const std::filesystem::path& path) {
        auto field = std::make_shared<traccc::magnetic_field>();
        traccc::io::read_magnetic_field(*field, path.string());
        return field;
      },
      "path"_a);

  auto f =
      py::class_<TracccMagneticField, std::shared_ptr<TracccMagneticField>>(
          m, "TracccMagneticField")
          .def(py::init([](const TracccMagneticField::Config& cfg,
                           Acts::Logging::Level level) {
                 return std::make_shared<TracccMagneticField>(
                     cfg, Acts::getDefaultLogger("TracccMagneticField", level));
               }),
               "config"_a, "level"_a = Acts::Logging::INFO)
          .def_property_readonly("config", &TracccMagneticField::config);

  auto c =
      py::class_<TracccMagneticField::Config>(f, "Config").def(py::init<>());
  ACTS_PYTHON_STRUCT(c, backend, field);
}

void addCopyTool(py::module_& m) {
  py::class_<TracccCopyTool, std::shared_ptr<TracccCopyTool>>(m,
                                                              "TracccCopyTool")
      .def(py::init<std::shared_ptr<const TracccBackend>>(), "backend"_a);
}

void addAlgorithms(py::module_& m) {
  // Host <-> device copies, one instance per traccc collection
  {
    using Alg = TracccHostToDeviceCopyAlg<traccc::edm::measurement_collection>;
    ACTS_PYTHON_DECLARE_ALGORITHM(Alg, m,
                                  "TracccHostToDeviceMeasurementCopyAlg",
                                  copyTool, input, output);
  }
  {
    using Alg = TracccDeviceToHostCopyAlg<traccc::edm::measurement_collection>;
    ACTS_PYTHON_DECLARE_ALGORITHM(Alg, m,
                                  "TracccDeviceToHostMeasurementCopyAlg",
                                  copyTool, input, output, inputMeasurements);
  }
  {
    using Alg = TracccDeviceToHostCopyAlg<traccc::edm::spacepoint_collection>;
    ACTS_PYTHON_DECLARE_ALGORITHM(Alg, m, "TracccDeviceToHostSpacepointCopyAlg",
                                  copyTool, input, output, inputMeasurements);
  }
  {
    using Alg = TracccDeviceToHostCopyAlg<traccc::edm::seed_collection>;
    ACTS_PYTHON_DECLARE_ALGORITHM(Alg, m, "TracccDeviceToHostSeedCopyAlg",
                                  copyTool, input, output, inputMeasurements);
  }
  {
    using Alg = TracccDeviceToHostCopyAlg<
        traccc::bound_track_parameters_collection_types>;
    ACTS_PYTHON_DECLARE_ALGORITHM(Alg, m,
                                  "TracccDeviceToHostTrackParametersCopyAlg",
                                  copyTool, input, output, inputMeasurements);
  }
  {
    using Alg = TracccDeviceToHostCopyAlg<TracccTrackContainer>;
    ACTS_PYTHON_DECLARE_ALGORITHM(Alg, m, "TracccDeviceToHostTrackCopyAlg",
                                  copyTool, input, output, inputMeasurements);
  }
}

}  // namespace

PYBIND11_MODULE(ActsExamplesPythonBindingsTraccc, traccc) {
  addBackends(traccc);
  addDetector(traccc);
  addMagneticField(traccc);
  addCopyTool(traccc);
  addAlgorithms(traccc);
}

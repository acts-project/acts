// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/Traccc/TracccBackend.hpp"
#include "ActsExamples/Traccc/TracccCopyAlgs.hpp"
#include "ActsExamples/Traccc/TracccCopyTool.hpp"
#include "ActsExamples/Traccc/TracccDetector.hpp"
#include "ActsExamples/Traccc/TracccEventData.hpp"
#include "ActsPython/Utilities/Helpers.hpp"
#include "ActsPython/Utilities/Macros.hpp"

#ifdef ACTS_EXAMPLES_TRACCC_CUDA
#include "ActsExamples/Traccc/Cuda/CudaTracccBackend.hpp"
#endif

#include <memory>
#include <string>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

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
                     detrayDetector, constantBField, bFieldFile);
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
  addCopyTool(traccc);
  addAlgorithms(traccc);
}

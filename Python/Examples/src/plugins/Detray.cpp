// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/Detray/DetrayPropagator.hpp"
#include "ActsExamples/Propagation/PropagatorInterface.hpp"
#include "ActsPlugins/Covfie/FieldConversion.hpp"
#include "ActsPlugins/Detray/DetrayConversionUtils.hpp"
#include "ActsPython/Utilities/Helpers.hpp"

// detray includes
#include <detray/core/detector.hpp>
#include <detray/detectors/odd_metadata.hpp>
#include <detray/io/frontend/detector_reader.hpp>
#include <detray/navigation/volume_graph.hpp>
#include <detray/propagator/line_stepper.hpp>
#include <detray/propagator/propagator.hpp>
#include <detray/propagator/rk_stepper.hpp>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <vecmem/memory/host_memory_resource.hpp>
#include <vecmem/memory/memory_resource.hpp>

namespace py = pybind11;
using namespace pybind11::literals;

using namespace Acts;
using namespace ActsExamples;
using namespace ActsPlugins;
using namespace ActsPython;

PYBIND11_MODULE(ActsExamplesPythonBindingsDetray, detray) {
  /// Propagators for ODD
  {
    using DetrayMetaDataODD = detray::odd_metadata<detray::array<float>>;

    using DetrayDetectorODD = detray::detector<DetrayMetaDataODD>;

    detray.def(
        "StraightLinePropagatorODD",
        [](std::shared_ptr<DetrayDetectorODD> detrayDetector,
           vecmem::memory_resource& memoryResource, bool sterile,
           Acts::Logging::Level logLevel,
           std::array<unsigned int, 2> searchWindow) {
          std::shared_ptr<PropagatorInterface> detrayProagator = nullptr;

          using DetrayLineStepper =
              detray::line_stepper<typename DetrayDetectorODD::algebra_type>;

          if (sterile) {
            using DetrayPropagator =
                DetraySterilePropagator<DetrayLineStepper, DetrayDetectorODD>;

            detrayProagator = std::make_shared<DetrayPropagator>(
                detrayDetector, memoryResource,
                Acts::getDefaultLogger("DetrayPropagator", logLevel),
                searchWindow);
          } else {
            using DetrayPropagator =
                DetrayPropagator<DetrayLineStepper, DetrayDetectorODD>;

            detrayProagator = std::make_shared<DetrayPropagator>(
                detrayDetector, memoryResource,
                Acts::getDefaultLogger("DetrayPropagator", logLevel),
                searchWindow);
          }
          return detrayProagator;
        },
        py::arg("detrayDetector"), py::arg("memoryResource"),
        py::arg("sterile"), py::arg("logLevel") = Acts::Logging::INFO,
        py::arg("searchWindow") = std::array<unsigned int, 2>{0u, 0u});
  }
}

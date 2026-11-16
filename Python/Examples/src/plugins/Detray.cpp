// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/Detray/DetrayPropagator.hpp"
#include "ActsExamples/Propagation/PropagatorInterface.hpp"

#include <detray/propagator/line_stepper.hpp>
#include <detray/python/detector.hpp>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <vecmem/memory/memory_resource.hpp>

namespace py = pybind11;
using namespace ActsExamples;

PYBIND11_MODULE(ActsExamplesPythonBindingsDetray, m) {
  py::module_::import("detray.core");
  m.def(
      "StraightLinePropagator",
      [](std::shared_ptr<detray::host_detector> detector,
         vecmem::memory_resource& resource, bool sterile,
         Acts::Logging::Level logLevel, std::array<unsigned int, 2> window) {
        return detray::python::visit(
            *detector,
            [&](const auto& concrete) -> std::shared_ptr<PropagatorInterface> {
              using Detector = std::remove_cvref_t<decltype(concrete)>;
              using Stepper =
                  detray::line_stepper<typename Detector::algebra_type>;
              auto owner = std::shared_ptr<const Detector>(detector, &concrete);
              auto logger =
                  Acts::getDefaultLogger("DetrayPropagator", logLevel);
              if (sterile) {
                return std::make_shared<
                    DetraySterilePropagator<Stepper, Detector>>(
                    owner, resource, std::move(logger), window);
              }
              return std::make_shared<DetrayPropagator<Stepper, Detector>>(
                  owner, resource, std::move(logger), window);
            });
      },
      py::arg("detrayDetector").none(false),
      py::arg("memoryResource").none(false), py::arg("sterile"),
      py::arg("logLevel") = Acts::Logging::INFO,
      py::arg("searchWindow") = std::array<unsigned int, 2>{0u, 0u},
      py::keep_alive<0, 1>(), py::keep_alive<0, 2>());
}

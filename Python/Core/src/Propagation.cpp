// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Propagator/Propagator.hpp"

#include "Acts/MagneticField/MagneticFieldProvider.hpp"
#include "Acts/Navigation/DetectorNavigator.hpp"
#include "Acts/Propagator/AtlasStepper.hpp"
#include "Acts/Propagator/EigenStepper.hpp"
#include "Acts/Propagator/Navigator.hpp"
#include "Acts/Propagator/StraightLineStepper.hpp"
#include "Acts/Propagator/SympyStepper.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "Acts/Utilities/TypeTag.hpp"
#include "ActsPython/Utilities/Context.hpp"
#include "ActsPython/Utilities/Macros.hpp"

#include <memory>
#include <stdexcept>
#include <utility>

#include <boost/core/demangle.hpp>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;
using namespace pybind11::literals;
using namespace Acts;

namespace ActsPython {
void addPropagation(Context& ctx) {
  auto& m = ctx.get("main");

  {
    using Config = Navigator::Config;
    auto nav =
        py::class_<Navigator, std::shared_ptr<Navigator>>(m, "Navigator")
            .def(py::init<>([](Config cfg,
                               Logging::Level level = Logging::INFO) {
                   return Navigator{cfg, getDefaultLogger("Navigator", level)};
                 }),
                 py::arg("cfg"), py::arg("level") = Logging::INFO);

    auto c = py::class_<Config>(nav, "Config").def(py::init<>());

    ACTS_PYTHON_STRUCT(c, resolveMaterial, resolvePassive, resolveSensitive,
                       trackingGeometry);
  }

  {
    using Config = Experimental::DetectorNavigator::Config;
    auto nav =
        py::class_<Experimental::DetectorNavigator,
                   std::shared_ptr<Experimental::DetectorNavigator>>(
            m, "DetectorNavigator")
            .def(py::init<>(
                     [](Config cfg, Logging::Level level = Logging::INFO) {
                       return Experimental::DetectorNavigator{
                           cfg, getDefaultLogger("DetectorNavigator", level)};
                     }),
                 py::arg("cfg"), py::arg("level") = Logging::INFO);

    auto c = py::class_<Config>(nav, "Config").def(py::init<>());

    ACTS_PYTHON_STRUCT(c, resolveMaterial, resolvePassive, resolveSensitive,
                       detector);
  }

  {
    auto stepper = py::class_<Acts::AtlasStepper>(m, "AtlasStepper");
    stepper.def(py::init<std::shared_ptr<const Acts::MagneticFieldProvider>>());
  }

  {
    auto stepper = py::class_<EigenStepper<>>(m, "EigenStepper");
    stepper.def(py::init<std::shared_ptr<const MagneticFieldProvider>>());
  }

  {
    auto stepper = py::class_<Acts::SympyStepper>(m, "SympyStepper");
    stepper.def(py::init<std::shared_ptr<const Acts::MagneticFieldProvider>>());
  }

  {
     auto stepper =
        py::class_<Acts::StraightLineStepper>(m, "StraightLineStepper");
    stepper.def(py::init<>());
  }
}
}  // namespace ActsPython
// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/EventData/BoundTrackParameters.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/MagneticField/MagneticFieldContext.hpp"
#include "Acts/Propagator/AtlasStepper.hpp"
#include "Acts/Propagator/EigenStepper.hpp"
#include "Acts/Propagator/Navigator.hpp"
#include "Acts/Propagator/Propagator.hpp"
#include "Acts/Propagator/PropagatorOptions.hpp"
#include "Acts/Propagator/PropagatorTraits.hpp"
#include "Acts/Propagator/StraightLineStepper.hpp"
#include "Acts/Propagator/SympyStepper.hpp"
#include "Acts/Propagator/VoidNavigator.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "ActsPython/Utilities/Helpers.hpp"
#include "ActsPython/Utilities/Macros.hpp"

#include <memory>
#include <sstream>
#include <stdexcept>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;
using namespace pybind11::literals;

using namespace Acts;

namespace ActsPython {

/// @brief Propagate to a target surface, throwing an exception if the
/// result is invalid.
BoundTrackParameters unwrapPropagateToSurface(
    Result<BoundTrackParameters> res) {
  if (!res.ok()) {
    std::stringstream ss;
    ss << "Propagation to surface failed with error: \"" << res.error() << "\"";
    throw std::runtime_error{ss.str()};
  }
  return *res;
}

template <typename stepper_t, typename navigator_t>
void addPropagator(py::module_& m, const std::string& prefix) {
  using propagator_t = Propagator<stepper_t, navigator_t>;
  auto p =
      py::class_<propagator_t>(m, (prefix + "Propagator").c_str())
          .def(py::init<>([=](stepper_t stepper, navigator_t navigator,
                              Logging::Level level = Logging::Level::INFO) {
                 return propagator_t{
                     std::move(stepper), std::move(navigator),
                     getDefaultLogger(prefix + "Propagator", level)};
               }),
               py::arg("stepper"), py::arg("navigator"),
               py::arg("level") = Logging::INFO);

  if constexpr (SupportsBoundParameters_v<stepper_t>) {
    p.def(
        "propagateToSurface",
        [](const propagator_t& self, const BoundTrackParameters& start,
           const Surface& target, const PropagatorPlainOptions& options) {
          return unwrapPropagateToSurface(
              self.propagateToSurface(start, target, options));
        },
        py::arg("start"), py::arg("target"), py::arg("options"));
  }
}

/// @brief Bind propagation related functions to the Python module
/// @param m The Python module to which the functions will be bound
void addPropagation(py::module_& m) {
  {
    py::class_<PropagatorPlainOptions>(m, "PropagatorPlainOptions")
        .def(py::init<const GeometryContext&, const MagneticFieldContext&>(),
             py::arg("geoContext"), py::arg("magFieldContext"))
        .def_readwrite("maxSteps", &PropagatorPlainOptions::maxSteps)
        .def_readwrite("pathLimit", &PropagatorPlainOptions::pathLimit)
        .def_readwrite("loopProtection",
                       &PropagatorPlainOptions::loopProtection)
        .def_readwrite("surfaceTolerance",
                       &PropagatorPlainOptions::surfaceTolerance);
  }

  {
    using Config = Navigator::Config;
    auto nav =
        py::class_<Navigator, std::shared_ptr<Navigator>>(m, "Navigator")
            .def(py::init<>([](const Config& cfg,
                               Logging::Level level = Logging::INFO) {
                   return Navigator{cfg, getDefaultLogger("Navigator", level)};
                 }),
                 py::arg("cfg"), py::arg("level") = Logging::INFO);

    auto c = py::class_<Config>(nav, "Config").def(py::init<>());

    ACTS_PYTHON_STRUCT(c, resolveMaterial, resolvePassive, resolveSensitive,
                       trackingGeometry);
  }

  {
    py::class_<VoidNavigator, std::shared_ptr<VoidNavigator>>(m,
                                                              "VoidNavigator")
        .def(py::init<>());
  }

  {
    auto stepper = py::class_<AtlasStepper>(m, "AtlasStepper");
    stepper.def(py::init<std::shared_ptr<const MagneticFieldProvider>>());
    addPropagator<AtlasStepper, Navigator>(m, "Atlas");
  }
  {
    auto stepper = py::class_<EigenStepper<>>(m, "EigenStepper");
    stepper.def(py::init<std::shared_ptr<const MagneticFieldProvider>>());
    addPropagator<EigenStepper<>, Navigator>(m, "Eigen");
    addPropagator<EigenStepper<>, VoidNavigator>(m, "EigenVoid");
  }
  {
    auto stepper = py::class_<StraightLineStepper>(m, "StraightLineStepper");
    stepper.def(py::init<>());
    addPropagator<StraightLineStepper, Navigator>(m, "StraightLine");
  }
  {
    auto stepper = py::class_<SympyStepper>(m, "SympyStepper");
    stepper.def(py::init<std::shared_ptr<const MagneticFieldProvider>>());
    addPropagator<SympyStepper, Navigator>(m, "Sympy");
  }
}

}  // namespace ActsPython

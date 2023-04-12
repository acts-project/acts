// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Navigation/NextNavigator.hpp"
#include "Acts/Plugins/Python/Utilities.hpp"
#include "Acts/Propagator/AtlasStepper.hpp"
#include "Acts/Propagator/EigenStepper.hpp"
#include "Acts/Propagator/Navigator.hpp"
#include "Acts/Propagator/Propagator.hpp"
#include "Acts/Propagator/StraightLineStepper.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/Propagation/PropagationAlgorithm.hpp"
#include "ActsExamples/Propagation/PropagatorInterface.hpp"

#include <memory>
#include <optional>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;

namespace {

template <typename stepper_t, typename navigator_t>
void addPropagator(py::module_& m, const std::string& prefix) {
  using propagator_t = Acts::Propagator<stepper_t, navigator_t>;
  py::class_<propagator_t>(m, (prefix + "Propagator").c_str())
      .def(py::init<>(
               [=](stepper_t stepper, navigator_t navigator,
                   Acts::Logging::Level level = Acts::Logging::Level::INFO) {
                 return propagator_t{
                     std::move(stepper), std::move(navigator),
                     Acts::getDefaultLogger(prefix + "Propagator", level)};
               }),
           py::arg("stepper"), py::arg("navigator"),
           py::arg("level") = Acts::Logging::INFO);

  using prop_if_t = ActsExamples::ConcretePropagator<propagator_t>;
  py::class_<prop_if_t, ActsExamples::PropagatorInterface,
             std::shared_ptr<prop_if_t>>(
      m, (prefix + "ConcretePropagator").c_str())
      .def(py::init<propagator_t>());
}

}  // namespace

namespace Acts::Python {
void addPropagation(Context& ctx) {
  auto [m, prop, mex] = ctx.get("main", "propagation", "examples");

  {
    using Config = Acts::Navigator::Config;
    auto nav =
        py::class_<Acts::Navigator, std::shared_ptr<Acts::Navigator>>(
            m, "Navigator")
            .def(py::init<>([](Config cfg,
                               Logging::Level level = Logging::INFO) {
                   return Navigator{cfg, getDefaultLogger("Navigator", level)};
                 }),
                 py::arg("cfg"), py::arg("level") = Logging::INFO);

    auto c = py::class_<Config>(nav, "Config").def(py::init<>());

    ACTS_PYTHON_STRUCT_BEGIN(c, Config);
    ACTS_PYTHON_MEMBER(resolveMaterial);
    ACTS_PYTHON_MEMBER(resolvePassive);
    ACTS_PYTHON_MEMBER(resolveSensitive);
    ACTS_PYTHON_MEMBER(trackingGeometry);
    ACTS_PYTHON_STRUCT_END();
  }

  {
    using Config = Acts::Experimental::NextNavigator::Config;
    auto nav = py::class_<Acts::Experimental::NextNavigator,
                          std::shared_ptr<Acts::Experimental::NextNavigator>>(
                   m, "NextNavigator")
                   .def(py::init<>([](Config cfg,
                                      Logging::Level level = Logging::INFO) {
                          return Acts::Experimental::NextNavigator{
                              cfg, getDefaultLogger("NextNavigator", level)};
                        }),
                        py::arg("cfg"), py::arg("level") = Logging::INFO);

    auto c = py::class_<Config>(nav, "Config").def(py::init<>());

    ACTS_PYTHON_STRUCT_BEGIN(c, Config);
    ACTS_PYTHON_MEMBER(resolveMaterial);
    ACTS_PYTHON_MEMBER(resolvePassive);
    ACTS_PYTHON_MEMBER(resolveSensitive);
    ACTS_PYTHON_MEMBER(detector);
    ACTS_PYTHON_STRUCT_END();
  }

  ACTS_PYTHON_DECLARE_ALGORITHM(
      ActsExamples::PropagationAlgorithm, mex, "PropagationAlgorithm",
      propagatorImpl, randomNumberSvc, mode, sterileLogger, debugOutput,
      energyLoss, multipleScattering, recordMaterialInteractions, ntests,
      d0Sigma, z0Sigma, phiSigma, thetaSigma, qpSigma, tSigma, phiRange,
      etaRange, ptRange, ptLoopers, maxStepSize, propagationStepCollection,
      propagationMaterialCollection, covarianceTransport, covariances,
      correlations);

  py::class_<ActsExamples::PropagatorInterface,
             std::shared_ptr<ActsExamples::PropagatorInterface>>(
      mex, "PropagatorInterface");

  {
    auto stepper = py::class_<Acts::EigenStepper<>>(m, "EigenStepper");
    stepper.def(
        // Add custom constructor lambda so that not specifying the overstep
        // limit takes the default from C++ EigenStepper
        py::init(
            [](std::shared_ptr<const Acts::MagneticFieldProvider> bField,
               std::optional<double> overstepLimit) -> Acts::EigenStepper<> {
              if (overstepLimit) {
                return {std::move(bField), overstepLimit.value()};
              } else {
                return {std::move(bField)};
              }
            }),
        py::arg("bField"), py::arg("overstepLimit") = std::nullopt);

    addPropagator<Acts::EigenStepper<>, Acts::Navigator>(prop, "Eigen");
  }

  {
    addPropagator<Acts::EigenStepper<>, Acts::Experimental::NextNavigator>(
        prop, "EigenNext");
  }

  {
    auto stepper = py::class_<Acts::AtlasStepper>(m, "AtlasStepper");
    stepper.def(py::init<std::shared_ptr<const Acts::MagneticFieldProvider>>());

    addPropagator<Acts::AtlasStepper, Acts::Navigator>(prop, "Atlas");
  }

  {
    auto stepper =
        py::class_<Acts::StraightLineStepper>(m, "StraightLineStepper");
    stepper.def(py::init<>());

    addPropagator<Acts::StraightLineStepper, Acts::Navigator>(prop,
                                                              "StraightLine");
  }
}

}  // namespace Acts::Python

// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

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
template <typename stepper_t>
auto addStepper(const std::string& prefix, py::module_& m, py::module_& prop) {
  auto stepper = py::class_<stepper_t>(m, (prefix + "Stepper").c_str());

  using propagator_t = Acts::Propagator<stepper_t, Acts::Navigator>;
  auto propagator =
      py::class_<propagator_t>(prop, (prefix + "Propagator").c_str())
          .def(py::init<>([=](stepper_t _stepper, Acts::Navigator navigator,
                              Acts::Logging::Level level =
                                  Acts::Logging::Level::INFO) {
                 return propagator_t{
                     std::move(_stepper), std::move(navigator),
                     Acts::getDefaultLogger(prefix + "Propagator", level)};
               }),
               py::arg("stepper"), py::arg("navigator"),
               py::arg("level") = Acts::Logging::INFO);

  using prop_if_t = ActsExamples::ConcretePropagator<propagator_t>;
  py::class_<prop_if_t, ActsExamples::PropagatorInterface,
             std::shared_ptr<prop_if_t>>(
      prop, (prefix + "ConcretePropagator").c_str())
      .def(py::init<propagator_t>());

  return std::pair{stepper, propagator};
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
    auto [stepper, propagator] =
        addStepper<Acts::EigenStepper<>>("Eigen", m, prop);
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
  }

  {
    auto [stepper, propagator] =
        addStepper<Acts::AtlasStepper>("Atlas", m, prop);
    stepper.def(py::init<std::shared_ptr<const Acts::MagneticFieldProvider>>());
  }

  {
    auto [stepper, propagator] =
        addStepper<Acts::StraightLineStepper>("StraightLine", m, prop);
    stepper.def(py::init<>());
  }
}

}  // namespace Acts::Python

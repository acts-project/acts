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
#include "ActsExamples/Propagation/PropagationAlgorithm.hpp"
#include "ActsExamples/Propagation/PropagatorInterface.hpp"

#include <memory>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;

namespace {
template <typename stepper_t>
void addStepper(const std::string& prefix, py::module_& m, py::module_& prop) {
  auto stepper = py::class_<stepper_t>(m, (prefix + "Stepper").c_str());

  if constexpr (std::is_same_v<stepper_t, Acts::StraightLineStepper>) {
    stepper.def(py::init<>());
  } else {
    stepper.def(py::init<std::shared_ptr<const Acts::MagneticFieldProvider>>());
  }

  using propagator_t = Acts::Propagator<stepper_t, Acts::Navigator>;
  py::class_<propagator_t>(prop, (prefix + "Propagator").c_str())
      .def(py::init<stepper_t, Acts::Navigator>());

  using prop_if_t = ActsExamples::ConcretePropagator<propagator_t>;
  py::class_<prop_if_t, ActsExamples::PropagatorInterface,
             std::shared_ptr<prop_if_t>>(
      prop, (prefix + "ConcretePropagator").c_str())
      .def(py::init<propagator_t>());
}

}  // namespace

namespace Acts::Python {
void addPropagation(Context& ctx) {
  auto [m, prop, mex] = ctx.get("main", "propagation", "examples");

  {
    using Config = Acts::Navigator::Config;
    auto nav = py::class_<Acts::Navigator, std::shared_ptr<Acts::Navigator>>(
                   m, "Navigator")
                   .def(py::init<Config>());

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

  addStepper<Acts::EigenStepper<>>("Eigen", m, prop);
  addStepper<Acts::AtlasStepper>("Atlas", m, prop);
  addStepper<Acts::StraightLineStepper>("StraightLine", m, prop);
}

}  // namespace Acts::Python

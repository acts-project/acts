// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Definitions/Direction.hpp"
#include "Acts/Propagator/AtlasStepper.hpp"
#include "Acts/Propagator/EigenStepper.hpp"
#include "Acts/Propagator/Navigator.hpp"
#include "Acts/Propagator/Propagator.hpp"
#include "Acts/Propagator/StraightLineStepper.hpp"
#include "Acts/Propagator/SympyStepper.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/Propagation/PropagationAlgorithm.hpp"
#include "ActsExamples/Propagation/PropagatorInterface.hpp"
#include "ActsPython/Utilities/Helpers.hpp"
#include "ActsPython/Utilities/Macros.hpp"

#include <algorithm>
#include <array>
#include <map>
#include <memory>
#include <optional>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace Acts {
class MagneticFieldProvider;
}  // namespace Acts

namespace py = pybind11;
using namespace Acts;
using namespace ActsExamples;

namespace {

template <typename stepper_t, typename navigator_t>
void addConcretePropagator(py::module_& m, const std::string& prefix) {
  using propagator_t = Propagator<stepper_t, navigator_t>;

  using prop_if_t = ConcretePropagator<propagator_t>;
  py::class_<prop_if_t, PropagatorInterface, std::shared_ptr<prop_if_t>>(
      m, (prefix + "ConcretePropagator").c_str())
      .def(py::init<propagator_t>());
}

}  // namespace

namespace ActsPython {
void addPropagation(py::module& mex) {
  ACTS_PYTHON_DECLARE_ALGORITHM(
      PropagationAlgorithm, mex, "PropagationAlgorithm", propagatorImpl,
      sterileLogger, debugOutput, energyLoss, multipleScattering,
      recordMaterialInteractions, ptLoopers, maxStepSize, covarianceTransport,
      inputTrackParameters, outputSummaryCollection, outputMaterialCollection);

  py::class_<PropagatorInterface, std::shared_ptr<PropagatorInterface>>(
      mex, "PropagatorInterface");

  // Eigen stepper based propagator
  { addConcretePropagator<EigenStepper<>, Navigator>(mex, "Eigen"); }

  // ATLAS stepper based propagator
  { addConcretePropagator<AtlasStepper, Navigator>(mex, "Atlas"); }

  // Sympy stepper based propagator
  { addConcretePropagator<SympyStepper, Navigator>(mex, "Sympy"); }

  // Straight line stepper
  {
    addConcretePropagator<StraightLineStepper, Navigator>(mex, "StraightLine");
  }
}

}  // namespace ActsPython

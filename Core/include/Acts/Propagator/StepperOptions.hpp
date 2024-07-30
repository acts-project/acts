// This file is part of the Acts project.
//
// Copyright (C) 2024 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Units.hpp"

#include <limits>

namespace Acts {

struct StepperPlainOptions {
  /// Tolerance for the error of the integration
  double stepTolerance = 1e-4;

  /// Cut-off value for the step size
  double stepSizeCutOff = 0.;

  /// Absolute maximum step size
  double maxStepSize = 10 * Acts::UnitConstants::m;

  /// Maximum number of Runge-Kutta steps for the stepper step call
  unsigned int maxRungeKuttaStepTrials = 10000;

  struct Dense {
    /// Toggle between mean and mode evaluation of energy loss
    bool meanEnergyLoss = true;

    /// Boolean flag for inclusion of d(dEds)d(q/p) into energy loss
    bool includeGradient = true;

    /// Cut-off value for the momentum in SI units
    double momentumCutOff = 0.;
  } dense;
};

}  // namespace Acts

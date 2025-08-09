// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Propagator/EigenStepper.hpp"
#include "Acts/Propagator/MultiStepperLoop.hpp"

namespace Acts {

/// Alias for the multi-stepper based on the Acts::EigenStepper
template <typename extension_t = EigenStepperDefaultExtension,
          typename reducer_t = MaxWeightReducerLoop>
using MultiEigenStepperLoop =
    MultiStepperLoop<EigenStepper<extension_t>, reducer_t>;

}  // namespace Acts

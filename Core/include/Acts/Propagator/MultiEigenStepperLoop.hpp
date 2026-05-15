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

/// @brief Type alias for loop reducer based on maximum momentum
/// @details Reduces multiple loop iterations based on momentum criteria
using MaxMomentumReducerLoop =
    detail::SingleComponentReducer<detail::MaxMomentumComponent>;

/// @brief Type alias for loop reducer based on maximum weight
/// @details Reduces multiple loop iterations based on weight criteria
using MaxWeightReducerLoop =
    detail::SingleComponentReducer<detail::MaxWeightComponent>;

/// @brief Stepper based on the EigenStepper, but handles Multi-Component Tracks
/// (e.g., for the GSF). Internally, this only manages a vector of
/// EigenStepper::States. This simplifies implementation, but has several
/// drawbacks:
/// * There are certain redundancies between the global State and the
/// component states
/// * The components do not share a single magnetic-field-cache
/// @tparam extension_t See EigenStepper for details
/// @tparam component_reducer_t How to map the multi-component state to a single
/// component
/// @tparam small_vector_size A size-hint how much memory should be allocated
/// by the small vector
template <typename extension_t = EigenStepperDefaultExtension,
          typename reducer_t = MaxWeightReducerLoop>
using MultiEigenStepperLoop =
    MultiStepperLoop<EigenStepper<extension_t>, reducer_t>;

}  // namespace Acts

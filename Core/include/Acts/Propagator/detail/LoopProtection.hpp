// This file is part of the Acts project.
//
// Copyright (C) 2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Utilities/Helpers.hpp"

namespace Acts {
namespace detail {

template <typename path_arborter_t>
struct LoopProtection {
  /// @brief Call to dress the options with a loop protection
  ///
  /// @tparam state_t State type of the Propagation call
  /// @tparam stepper_t Stepper type of the Propagator setup
  ///
  /// @param [in,out] state State object provided for the call
  /// @param [in] stepper Stepper used
  template <typename propagator_state_t, typename stepper_t>
  void operator()(propagator_state_t& state, const stepper_t& stepper) const {
    const auto& logger = state.options.logger;
    // Estimate the loop protection limit
    if (state.options.loopProtection) {
      // Get the field at the start position
      Vector3 field =
          stepper.getField(state.stepping, stepper.position(state.stepping));
      const double B = field.norm();
      if (B != 0.) {
        // Transverse component at start is taken for the loop protection
        const double p = stepper.momentum(state.stepping);
        // Calculate the full helix path
        const double helixPath = state.stepping.navDir * 2 * M_PI * p / B;
        // And set it as the loop limit if it overwrites the internal limit
        auto& pathAborter =
            state.options.abortList.template get<path_arborter_t>();
        double loopLimit = state.options.loopFraction * helixPath;
        double pathLimit = pathAborter.internalLimit;
        if (std::abs(loopLimit) < std::abs(pathLimit)) {
          pathAborter.internalLimit = loopLimit;

          ACTS_VERBOSE("Path aborter limit set to "
                       << loopLimit << " (full helix =  " << helixPath << ")");
        }
      }
    }
  }
};

}  // namespace detail
}  // namespace Acts

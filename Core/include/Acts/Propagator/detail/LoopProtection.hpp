// This file is part of the Acts project.
//
// Copyright (C) 2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Utilities/Definitions.hpp"
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
    // Estimate the loop protection limit
    if (state.options.loopProtection) {
      // Get the field at the start position
      Vector3D field =
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
        if (loopLimit * loopLimit < pathLimit * pathLimit) {
          pathAborter.internalLimit = loopLimit;
          debugLog(state, [&] {
            std::stringstream dstream;
            dstream << "Path aborter limit set to ";
            dstream << loopLimit << " (full helix =  " << helixPath << ")";
            return dstream.str();
          });
        }
      }
    }
  }

  /// The private loop protection debug logging
  ///
  /// It needs to be fed by a lambda function that returns a string,
  /// that guarantees that the lambda is only called in the
  /// state.options.debug == true
  /// case in order not to spend time when not needed.
  ///
  /// @tparam propagator_state_t Type of the propagator state
  ///
  /// @param[in,out] state the propagator state for the debug flag,
  ///      prefix and length
  /// @param logAction is a callable function that returns a streamable object
  template <typename propagator_state_t>
  void debugLog(propagator_state_t& state,
                const std::function<std::string()>& logAction) const {
    if (state.options.debug) {
      std::vector<std::string> lines;
      std::string input = logAction();
      boost::split(lines, input, boost::is_any_of("\n"));
      for (const auto& line : lines) {
        std::stringstream dstream;
        dstream << '\n';
        dstream << " ∞ " << std::setw(state.options.debugPfxWidth)
                << " loop protection "
                << " | ";
        dstream << std::setw(state.options.debugMsgWidth) << line << '\n';
        state.options.debugString += dstream.str();
      }
    }
  }
};

}  // namespace detail
}  // namespace Acts

// This file is part of the Acts project.
//
// Copyright (C) 2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Utilities/Units.hpp"

namespace Acts {

namespace detail {

  template <typename path_arborter_t>
  struct LoopProtection
  {

    /// @brief Call to dress the options with a loop protection
    ///
    /// @tparam state_t State type of the Propagation call
    /// @tparam stepper_t Stepper type of the Propagator setup
    ///
    /// @param [in,out] state State object provided for the call
    /// @param [in] stepper Stepper used
    template <typename state_t, typename stepper_t>
    void
    operator()(state_t& state, const stepper_t& stepper) const
    {
      // Estimate the loop protection limit
      if (state.options.loopProtection) {
        // Get the field at the start position
        Vector3D field = stepper.getField(state.stepping,
                                          stepper.position(state.stepping));
        // Momentum in SI units and B field
        const double p
            = units::Nat2SI<units::MOMENTUM>(stepper.momentum(state.stepping));
        const double B         = field.norm();
        const double helixPath = 2 * M_PI * p / B;
        // now set the new loop limit
        auto& pathAborter
            = state.options.abortList.template get<path_arborter_t>();
        pathAborter.internalLimit = state.options.loopFraction * helixPath;
      }
    }
  };

}  // namespace detail
}  // namespace Acts

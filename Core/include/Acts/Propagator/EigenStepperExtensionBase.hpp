// This file is part of the Acts project.
//
// Copyright (C) 2024 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"

#include <array>

namespace Acts {

template <typename scalar_t, typename derived_t>
struct EigenStepperExtensionBase {
  using Scalar = scalar_t;
  /// @brief Vector3 replacement for the custom scalar type
  using ThisVector3 = Eigen::Matrix<Scalar, 3, 1>;

  /// @brief This functions broadcasts the call for evaluating k1. It collects
  /// all arguments and extensions, test their validity for the evaluation and
  /// passes them forward for evaluation and returns a boolean as indicator if
  /// the evaluation is valid.
  template <typename propagator_state_t, typename stepper_t,
            typename navigator_t>
  bool k1(const propagator_state_t& state, const stepper_t& stepper,
          const navigator_t& navigator, Vector3& knew, const Vector3& bField,
          std::array<double, 4>& kQoP) {
    return static_cast<derived_t*>(this)->k(state, stepper, navigator, knew,
                                            bField, kQoP, 0);
  }

  /// @brief This functions broadcasts the call for evaluating k2. It collects
  /// all arguments and extensions and passes them forward for evaluation and
  /// returns a boolean as indicator if the evaluation is valid.
  template <typename propagator_state_t, typename stepper_t,
            typename navigator_t>
  bool k2(const propagator_state_t& state, const stepper_t& stepper,
          const navigator_t& navigator, Vector3& knew, const Vector3& bField,
          std::array<double, 4>& kQoP, const double h, const Vector3& kprev) {
    return static_cast<derived_t*>(this)->k(state, stepper, navigator, knew,
                                            bField, kQoP, 1, h, kprev);
  }

  /// @brief This functions broadcasts the call for evaluating k3. It collects
  /// all arguments and extensions and passes them forward for evaluation and
  /// returns a boolean as indicator if the evaluation is valid.
  template <typename propagator_state_t, typename stepper_t,
            typename navigator_t>
  bool k3(const propagator_state_t& state, const stepper_t& stepper,
          const navigator_t& navigator, Vector3& knew, const Vector3& bField,
          std::array<double, 4>& kQoP, const double h, const Vector3& kprev) {
    return static_cast<derived_t*>(this)->k(state, stepper, navigator, knew,
                                            bField, kQoP, 2, h, kprev);
  }

  /// @brief This functions broadcasts the call for evaluating k4. It collects
  /// all arguments and extensions and passes them forward for evaluation and
  /// returns a boolean as indicator if the evaluation is valid.
  template <typename propagator_state_t, typename stepper_t,
            typename navigator_t>
  bool k4(const propagator_state_t& state, const stepper_t& stepper,
          const navigator_t& navigator, Vector3& knew, const Vector3& bField,
          std::array<double, 4>& kQoP, const double h, const Vector3& kprev) {
    return static_cast<derived_t*>(this)->k(state, stepper, navigator, knew,
                                            bField, kQoP, 3, h, kprev);
  }
};

}  // namespace Acts

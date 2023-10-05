// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/Propagator/detail/Auctioneer.hpp"
#include "Acts/Utilities/detail/Extendable.hpp"
#include "Acts/Utilities/detail/MPL/all_of.hpp"
#include "Acts/Utilities/detail/MPL/has_duplicates.hpp"

#include <array>

namespace Acts {

/// @brief Container of extensions used in the stepper of the propagation. This
/// struct allows a broadcast of function calls for each element in the list.
/// The broadcasts occur for a certain function at each step in a specific
/// order.
/// The first function is an evaluater if an extension is or how many extensions
/// are applicable for an upcoming step.
/// The next functions called are the evaluations of the k_1 - k_4 or the RKN4
/// integration.
/// The last function call in a step is the finalize() method. This method is an
/// overloaded function (optionally propagates the covariance).
/// Each method has the possibility to break the evaluation of a given step if
/// an extension reports that something went wrong (e.g. a particle lost too
/// much momentum during the step)
/// @tparam extensions Types of the extensions
template <typename... extensions>
struct StepperExtensionList : private detail::Extendable<extensions...> {
 private:
  // Checkout for duplicates in the extensions
  static_assert(not detail::has_duplicates_v<extensions...>,
                "same extension type specified several times");

  static constexpr unsigned int nExtensions = sizeof...(extensions);

  static_assert(nExtensions != 0, "no extension type specified");

  // Access to all extensions
  using detail::Extendable<extensions...>::tuple;

  // Vector of valid extensions for a step
  std::array<bool, nExtensions> validExtensions{};

 public:
  // Access to an extension
  using detail::Extendable<extensions...>::get;

  /// @brief Evaluation function to set valid extensions for an upcoming
  /// integration step
  ///
  /// @tparam propagator_state_t Type of the state of the propagator
  /// @tparam stepper_t Type of the stepper
  /// @tparam navigtor_t Type of the navigator
  ///
  /// @param [in] state State of the propagator
  /// @param [in] stepper Stepper of the propagation
  /// @param [in] navigator Navigator of the propagation
  template <typename propagator_state_t, typename stepper_t,
            typename navigtor_t>
  bool validExtensionForStep(const propagator_state_t& state,
                             const stepper_t& stepper,
                             const navigtor_t& navigator) {
    const auto bids = std::apply(
        [&](const auto&... ext) {
          return std::array<int, nExtensions>{
              ext.bid(state, stepper, navigator)...};
        },
        tuple());

    validExtensions = state.stepping.auctioneer(std::move(bids));

    return (std::find(validExtensions.begin(), validExtensions.end(), true) !=
            validExtensions.end());
  }

  /// @brief This functions broadcasts the call for evaluating a generic k. It
  /// collects all arguments and extensions, test their validity for the
  /// evaluation and passes them forward for evaluation and returns a boolean as
  /// indicator if the evaluation is valid.
  template <typename propagator_state_t, typename stepper_t,
            typename navigator_t>
  bool k(const propagator_state_t& state, const stepper_t& stepper,
         const navigator_t& navigator, Vector3& knew, const Vector3& bField,
         std::array<double, 4>& kQoP, const int i, const double h = 0.,
         const Vector3& kprev = Vector3::Zero()) {
    // TODO replace with integer-templated lambda with C++20
    auto impl = [&, i, h](auto intType, auto& implRef) {
      constexpr int N = decltype(intType)::value;

      if constexpr (N == 0) {
        return true;
      } else {
        // If element is invalid: continue
        if (!std::get<N - 1>(validExtensions)) {
          return implRef(std::integral_constant<int, N - 1>{}, implRef);
        }
        // Continue as long as evaluations are 'true'
        if (std::get<N - 1>(this->tuple())
                .template k(state, stepper, navigator, knew, bField, kQoP, i, h,
                            kprev)) {
          return implRef(std::integral_constant<int, N - 1>{}, implRef);
        } else {
          // Break at false
          return false;
        }
      }
    };

    return impl(std::integral_constant<int, nExtensions>{}, impl);
  }

  /// @brief This functions broadcasts the call for evaluating k1. It collects
  /// all arguments and extensions, test their validity for the evaluation and
  /// passes them forward for evaluation and returns a boolean as indicator if
  /// the evaluation is valid.
  template <typename propagator_state_t, typename stepper_t,
            typename navigator_t>
  bool k1(const propagator_state_t& state, const stepper_t& stepper,
          const navigator_t& navigator, Vector3& knew, const Vector3& bField,
          std::array<double, 4>& kQoP) {
    return k(state, stepper, navigator, knew, bField, kQoP, 0);
  }

  /// @brief This functions broadcasts the call for evaluating k2. It collects
  /// all arguments and extensions and passes them forward for evaluation and
  /// returns a boolean as indicator if the evaluation is valid.
  template <typename propagator_state_t, typename stepper_t,
            typename navigator_t>
  bool k2(const propagator_state_t& state, const stepper_t& stepper,
          const navigator_t& navigator, Vector3& knew, const Vector3& bField,
          std::array<double, 4>& kQoP, const double h, const Vector3& kprev) {
    return k(state, stepper, navigator, knew, bField, kQoP, 1, h, kprev);
  }

  /// @brief This functions broadcasts the call for evaluating k3. It collects
  /// all arguments and extensions and passes them forward for evaluation and
  /// returns a boolean as indicator if the evaluation is valid.
  template <typename propagator_state_t, typename stepper_t,
            typename navigator_t>
  bool k3(const propagator_state_t& state, const stepper_t& stepper,
          const navigator_t& navigator, Vector3& knew, const Vector3& bField,
          std::array<double, 4>& kQoP, const double h, const Vector3& kprev) {
    return k(state, stepper, navigator, knew, bField, kQoP, 2, h, kprev);
  }

  /// @brief This functions broadcasts the call for evaluating k4. It collects
  /// all arguments and extensions and passes them forward for evaluation and
  /// returns a boolean as indicator if the evaluation is valid.
  template <typename propagator_state_t, typename stepper_t,
            typename navigator_t>
  bool k4(const propagator_state_t& state, const stepper_t& stepper,
          const navigator_t& navigator, Vector3& knew, const Vector3& bField,
          std::array<double, 4>& kQoP, const double h, const Vector3& kprev) {
    return k(state, stepper, navigator, knew, bField, kQoP, 3, h, kprev);
  }

  /// @brief This functions broadcasts the call of the method finalize(). It
  /// collects all extensions and arguments and passes them forward for
  /// evaluation and returns a boolean.
  template <typename propagator_state_t, typename stepper_t,
            typename navigator_t>
  bool finalize(propagator_state_t& state, const stepper_t& stepper,
                const navigator_t& navigator, const double h, FreeMatrix& D) {
    // TODO replace with integer-templated lambda with C++20
    auto impl = [&, h](auto intType, auto& implRef) {
      constexpr int N = decltype(intType)::value;

      if constexpr (N == 0) {
        return true;
      } else {
        // If element is invalid: continue
        if (!std::get<N - 1>(validExtensions)) {
          return implRef(std::integral_constant<int, N - 1>{}, implRef);
        }
        // Continue as long as evaluations are 'true'
        if (std::get<N - 1>(this->tuple())
                .finalize(state, stepper, navigator, h, D)) {
          return implRef(std::integral_constant<int, N - 1>{}, implRef);
        } else {
          // Break at false
          return false;
        }
      }
    };

    return impl(std::integral_constant<int, nExtensions>{}, impl);
  }

  /// @brief This functions broadcasts the call of the method finalize(). It
  /// collects all extensions and arguments and passes them forward for
  /// evaluation and returns a boolean.
  template <typename propagator_state_t, typename stepper_t,
            typename navigator_t>
  bool finalize(propagator_state_t& state, const stepper_t& stepper,
                const navigator_t& navigator, const double h) {
    // TODO replace with integer-templated lambda with C++20
    auto impl = [&, h](auto intType, auto& implRef) {
      constexpr int N = decltype(intType)::value;

      if constexpr (N == 0) {
        return true;
      } else {
        // If element is invalid: continue
        if (!std::get<N - 1>(validExtensions)) {
          return implRef(std::integral_constant<int, N - 1>{}, implRef);
        }

        // Continue as long as evaluations are 'true'
        if (std::get<N - 1>(this->tuple())
                .finalize(state, stepper, navigator, h)) {
          return implRef(std::integral_constant<int, N - 1>{}, implRef);
        } else {
          // Break at false
          return false;
        }
      }
    };

    return impl(std::integral_constant<int, nExtensions>{}, impl);
  }
};

}  // namespace Acts

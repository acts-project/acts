// This file is part of the Acts project.
//
// Copyright (C) 2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Propagator/detail/Auctioneer.hpp"
#include "Acts/Propagator/detail/stepper_extension_list_implementation.hpp"
#include "Acts/Utilities/detail/Extendable.hpp"
#include "Acts/Utilities/detail/MPL/all_of.hpp"
#include "Acts/Utilities/detail/MPL/has_duplicates.hpp"

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
struct StepperExtensionList : private detail::Extendable<extensions...>
{
private:
  // Checkout for duplicates in the extensions
  static_assert(not detail::has_duplicates_v<extensions...>,
                "same extension type specified several times");

  static constexpr unsigned int nExtensions = sizeof...(extensions);

  static_assert(nExtensions != 0, "no extension type specified");

  // Access to all extensions
  using detail::Extendable<extensions...>::tuple;

  using impl = detail::stepper_extension_list_impl<nExtensions>;

  // Vector of valid extensions for a step
  std::array<bool, nExtensions> validExtensions;

public:
  // Access to an extension
  using detail::Extendable<extensions...>::get;

  /// @brief Evaluation function to set valid extensions for an upcoming
  /// integration step
  ///
  /// @tparam propagator_state_t Type of the state of the propagator
  /// @tparam stepper_t Type of the stepper
  /// @param [in] state State of the propagator
  /// @param [in] stepper Stepper of the propagation
  template <typename propagator_state_t, typename stepper_t>
  bool
  validExtensionForStep(const propagator_state_t& state,
                        const stepper_t&          stepper)
  {
    std::array<int, nExtensions> bids;
    // Ask all extensions for a boolean statement of their validity
    impl::bid(tuple(), state, stepper, bids);
    // Post-process the vector in an auctioneer
    validExtensions = state.stepping.auctioneer(std::move(bids));

    return (std::find(validExtensions.begin(), validExtensions.end(), true)
            != validExtensions.end());
  }

  /// @brief This functions broadcasts the call for evaluating k1. It collects
  /// all arguments and extensions, test their validity for the evaluation and
  /// passes them forward for evaluation and returns a boolean as indicator if
  /// the evaluation is valid.
  template <typename propagator_state_t, typename stepper_t>
  bool
  k1(const propagator_state_t& state,
     const stepper_t&          stepper,
     Vector3D&                 knew,
     const Vector3D&           bField)
  {
    return impl::k(tuple(), state, stepper, knew, bField, validExtensions);
  }

  /// @brief This functions broadcasts the call for evaluating k2. It collects
  /// all arguments and extensions and passes them forward for evaluation and
  /// returns a boolean as indicator if the evaluation is valid.
  template <typename propagator_state_t, typename stepper_t>
  bool
  k2(const propagator_state_t& state,
     const stepper_t&          stepper,
     Vector3D&                 knew,
     const Vector3D&           bField,
     const double              h,
     const Vector3D&           kprev)
  {
    return impl::k(
        tuple(), state, stepper, knew, bField, validExtensions, 1, h, kprev);
  }

  /// @brief This functions broadcasts the call for evaluating k3. It collects
  /// all arguments and extensions and passes them forward for evaluation and
  /// returns a boolean as indicator if the evaluation is valid.
  template <typename propagator_state_t, typename stepper_t>
  bool
  k3(const propagator_state_t& state,
     const stepper_t&          stepper,
     Vector3D&                 knew,
     const Vector3D&           bField,
     const double              h,
     const Vector3D&           kprev)
  {
    return impl::k(
        tuple(), state, stepper, knew, bField, validExtensions, 2, h, kprev);
  }

  /// @brief This functions broadcasts the call for evaluating k4. It collects
  /// all arguments and extensions and passes them forward for evaluation and
  /// returns a boolean as indicator if the evaluation is valid.
  template <typename propagator_state_t, typename stepper_t>
  bool
  k4(const propagator_state_t& state,
     const stepper_t&          stepper,
     Vector3D&                 knew,
     const Vector3D&           bField,
     const double              h,
     const Vector3D&           kprev)
  {
    return impl::k(
        tuple(), state, stepper, knew, bField, validExtensions, 3, h, kprev);
  }

  /// @brief This functions broadcasts the call of the method finalize(). It
  /// collects all extensions and arguments and passes them forward for
  /// evaluation and returns a boolean.
  template <typename propagator_state_t, typename stepper_t>
  bool
  finalize(propagator_state_t& state,
           const stepper_t&    stepper,
           const double        h,
           ActsMatrixD<7, 7>& D)
  {
    return impl::finalize(tuple(), state, stepper, h, D, validExtensions);
  }

  /// @brief This functions broadcasts the call of the method finalize(). It
  /// collects all extensions and arguments and passes them forward for
  /// evaluation and returns a boolean.
  template <typename propagator_state_t, typename stepper_t>
  bool
  finalize(propagator_state_t& state, const stepper_t& stepper, const double h)
  {
    return impl::finalize(tuple(), state, stepper, h, validExtensions);
  }
};

}  // namespace Acts

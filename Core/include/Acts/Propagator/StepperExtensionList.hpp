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

/// @brief Container of extensions used in the stepper of the propagation
/// @tparam extensions Types of the extensions
template <typename... extensions>
struct StepperExtensionList : private detail::Extendable<extensions...>
{
private:
  // Checkout for duplicates in the extensions
  static_assert(not detail::has_duplicates_v<extensions...>,
                "same extension type specified several times");

  // Access to all extensions
  using detail::Extendable<extensions...>::tuple;

  using impl = detail::stepper_extension_list_impl<extensions...>;

  // Vector of valid extensions for a step
  std::vector<bool> validExtensions;

  /// @brief Evaluation function to set valid extensions for an upcoming
  /// integration step
  ///
  /// @tparam stepper_state_t Type of the state of the stepper
  /// @param [in] state State of the stepper
  template <typename stepper_state_t>
  void
  validExtensionForStep(const stepper_state_t& state)
  {
    std::vector<bool> validExtensionCandidates;
    // Ask all extensions for a boolean statement of their validity
    impl::validExtensionForStep(tuple(), state, validExtensionCandidates);
    // Post-process the vector in an auctioneer
    validExtensions = state.auctioneer(std::move(validExtensionCandidates));
  }

public:
  // Access to an extension
  using detail::Extendable<extensions...>::get;

  /// @brief This functions broadcasts the call of the method k(). It collects
  /// all extensions and arguments and passes them forward for evaluation and
  /// returns a boolean.
  template <typename stepper_state_t>
  bool
  k(const stepper_state_t& state,
    Vector3D&              knew,
    const Vector3D&        bField,
    const int              i     = 0,
    const double           h     = 0,
    const Vector3D&        kprev = Vector3D())
  {
    if (i == 0) {
      validExtensionForStep(state);
    }
    return impl::k(
        tuple(), state, knew, bField, validExtensions.cbegin(), i, h, kprev);
  }

  /// @brief This functions broadcasts the call of the method finalize(). It
  /// collects all extensions and arguments and passes them forward for
  /// evaluation and returns a boolean.
  template <typename stepper_state_t, typename stepper_data_t>
  bool
  finalize(stepper_state_t&      state,
           const double          h,
           const stepper_data_t& data,
           ActsMatrixD<7, 7>& D)
  {
    return impl::finalize(tuple(), state, h, data, D, validExtensions.cbegin());
  }

  /// @brief This functions broadcasts the call of the method finalize(). It
  /// collects all extensions and arguments and passes them forward for
  /// evaluation and returns a boolean.
  template <typename stepper_state_t>
  bool
  finalize(stepper_state_t& state, const double h)
  {
    return impl::finalize(tuple(), state, h, validExtensions.cbegin());
  }
};

}  // namespace Acts

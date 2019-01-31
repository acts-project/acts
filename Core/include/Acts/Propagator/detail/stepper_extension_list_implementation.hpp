// This file is part of the Acts project.
//
// Copyright (C) 2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once
#include "Acts/Utilities/Definitions.hpp"

/// @note The used variables in this file are:
///
/// @tparam N Size of the tuple @p obs_tuple
/// @tparam propagator_state_t Type of the state of the propagation
/// @tparam stepper_t Type of the stepper
/// @tparam T Types of the extensions
///
/// @param obs_tuple Extendable tuple of extensions
/// @param state State of the propagation
/// @param stepper Stepper of the propagation
/// @param bids List that collects bids about validity of an extension for an
/// upcoming step
/// @param validExtensions List of valid extensions that will be used
/// @param knew k_1 - k_4 that is about to be evaluated in the upcoming call
/// @param bField B-field at a certain point
/// @param i Defines the calculation of knew=k_{i+1}
/// @param h Distance to the starting point k_1
/// @param kprev k_i that is used for the evaluation for knew
/// @param D Transport matrix of the jacobian

namespace Acts {
namespace detail {
  /// The dummy list call implementation
  template <unsigned int N>
  struct stepper_extension_list_impl;

  /// The extension list call implementation
  /// - it calls 'validExtensionsForStep()' on the current entry of the tuple
  /// - it stores the result in @p validExtensions
  /// - then broadcasts the extension call to the remaining tuple
  template <unsigned int N>
  struct stepper_extension_list_impl
  {
    template <typename propagator_state_t, typename stepper_t, typename... T>
    static void
    bid(const std::tuple<T...>&   obs_tuple,
        const propagator_state_t& state,
        const stepper_t&          stepper,
        std::array<int, sizeof...(T)>& bids)
    {
      std::get<N - 1>(bids) = std::get<N - 1>(obs_tuple).bid(state, stepper);
      stepper_extension_list_impl<N - 1>::bid(obs_tuple, state, stepper, bids);
    }

    /// The extension list call implementation
    /// - it calls 'k()' on the current entry of the tuple
    /// - then broadcasts the extension call to the remaining tuple
    template <typename propagator_state_t, typename stepper_t, typename... T>
    static bool
    k(std::tuple<T...>&         obs_tuple,
      const propagator_state_t& state,
      const stepper_t&          stepper,
      Vector3D&                 knew,
      const Vector3D&           bField,
      const std::array<bool, sizeof...(T)>& validExtensions,
      const int       i     = 0,
      const double    h     = 0,
      const Vector3D& kprev = Vector3D())
    {
      // If element is invalid: continue
      if (!std::get<N - 1>(validExtensions)) {
        return stepper_extension_list_impl<N - 1>::k(obs_tuple,
                                                     state,
                                                     stepper,
                                                     knew,
                                                     bField,
                                                     validExtensions,
                                                     i,
                                                     h,
                                                     kprev);
      }

      // Continue as long as evaluations are 'true'
      if (std::get<N - 1>(obs_tuple).k(
              state, stepper, knew, bField, i, h, kprev)) {
        return stepper_extension_list_impl<N - 1>::k(obs_tuple,
                                                     state,
                                                     stepper,
                                                     knew,
                                                     bField,
                                                     validExtensions,
                                                     i,
                                                     h,
                                                     kprev);
      } else {
        // Break at false
        return false;
      }
    }

    /// The extension list call implementation
    /// - it calls 'finalize()' on the current entry of the tuple
    /// - then broadcasts the extension call to the remaining tuple
    template <typename propagator_state_t, typename stepper_t, typename... T>
    static bool
    finalize(const std::tuple<T...>& obs_tuple,
             propagator_state_t&     state,
             const stepper_t&        stepper,
             const double            h,
             ActsMatrixD<7, 7>&                    D,
             const std::array<bool, sizeof...(T)>& validExtensions)
    {
      // If element is invalid: continue
      if (!std::get<N - 1>(validExtensions)) {
        return stepper_extension_list_impl<N - 1>::finalize(
            obs_tuple, state, stepper, h, D, validExtensions);
      }

      // Continue as long as evaluations are 'true'
      if (std::get<N - 1>(obs_tuple).finalize(state, stepper, h, D)) {
        return stepper_extension_list_impl<N - 1>::finalize(
            obs_tuple, state, stepper, h, D, validExtensions);
      } else {
        // Break at false
        return false;
      }
    }

    /// The extension list call implementation
    /// - it calls 'finalize()' on the current entry of the tuple
    /// - then broadcasts the extension call to the remaining tuple
    template <typename propagator_state_t, typename stepper_t, typename... T>
    static bool
    finalize(const std::tuple<T...>& obs_tuple,
             propagator_state_t&     state,
             const stepper_t&        stepper,
             const double            h,
             const std::array<bool, sizeof...(T)>& validExtensions)
    {
      // If element is invalid: continue
      if (!std::get<N - 1>(validExtensions)) {
        return stepper_extension_list_impl<N - 1>::finalize(
            obs_tuple, state, stepper, h, validExtensions);
      }

      // Continue as long as evaluations are 'true'
      if (std::get<N - 1>(obs_tuple).finalize(state, stepper, h)) {
        return stepper_extension_list_impl<N - 1>::finalize(
            obs_tuple, state, stepper, h, validExtensions);
      } else {
        // Break at false
        return false;
      }
    }
  };

  /// Template specialized list call implementation
  template <>
  struct stepper_extension_list_impl<0u>
  {
    /// The empty extension list call implementation
    template <typename propagator_state_t, typename stepper_t, typename... T>
    static void
    bid(const std::tuple<T...>& /*unused*/,
        const propagator_state_t& /*unused*/,
        const stepper_t& /*unused*/,
        std::array<int, sizeof...(T)>& /*unused*/)
    {
    }

    /// The empty extension list call implementation
    template <typename propagator_state_t, typename stepper_t, typename... T>
    static bool
    k(std::tuple<T...>& /*unused*/,
      const propagator_state_t& /*unused*/,
      const stepper_t& /*unused*/,
      Vector3D& /*unused*/,
      const Vector3D& /*unused*/,
      const std::array<bool, sizeof...(T)>& /*unused*/,
      const int /*unused*/,
      const double /*unused*/,
      const Vector3D& /*unused*/)
    {
      return true;
    }

    /// The empty extension list call implementation
    template <typename propagator_state_t, typename stepper_t, typename... T>
    static bool
    finalize(const std::tuple<T...>& /*unused*/,
             propagator_state_t& /*unused*/,
             const stepper_t& /*unused*/,
             const double /*unused*/,
             ActsMatrixD<7, 7>& /*unused*/,
             const std::array<bool, sizeof...(T)>& /*unused*/)
    {
      return true;
    }

    /// The empty extension list call implementation
    template <typename propagator_state_t, typename stepper_t, typename... T>
    static bool
    finalize(const std::tuple<T...>& /*unused*/,
             propagator_state_t& /*unused*/,
             const stepper_t& /*unused*/,
             const double /*unused*/,
             const std::array<bool, sizeof...(T)>& /*unused*/)
    {
      return true;
    }
  };

}  // namespace detail
}  // namespace Acts

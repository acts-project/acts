// This file is part of the Acts project.
//
// Copyright (C) 2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once
#include "Acts/Utilities/Definitions.hpp"

namespace Acts {
namespace detail {
  /// The dummy list call implementation
  template <typename... extensions>
  struct extension_list_impl;

  /// The extension list call implementation
  /// - it calls 'k()' on the current entry of the tuple
  /// - then broadcasts the extension call to the remaining tuple
  template <typename first, typename... others>
  struct extension_list_impl<first, others...>
  {
    template <typename T, typename stepper_state_t>
    static bool
    k(T&                     obs_tuple,
      const stepper_state_t& state,
      Vector3D&              knew,
      const Vector3D&        bField,
      const int              i     = 0,
      const double           h     = 0,
      const Vector3D&        kprev = Vector3D())
    {
      auto& this_extension = std::get<first>(obs_tuple);
      // Continue as long as evaluations are 'true'
      if (this_extension.k(state, knew, bField, i, h, kprev))
        return extension_list_impl<others...>::k(
            obs_tuple, state, knew, bField, i, h, kprev);
      else
        // Break at false
        return false;
    }

    /// The extension list call implementation
    /// - it calls 'finalize()' on the current entry of the tuple
    /// - then broadcasts the extension call to the remaining tuple
    template <typename T, typename stepper_state_t, typename stepper_data_t>
    static bool
    finalize(const T&              obs_tuple,
             stepper_state_t&      state,
             const double          h,
             const stepper_data_t& data,
             ActsMatrixD<7, 7>& D)
    {
      auto& this_extension = std::get<first>(obs_tuple);
      // Continue as long as evaluations are 'true'
      if (this_extension.finalize(state, h, data, D))
        return extension_list_impl<others...>::finalize(
            obs_tuple, state, h, data, D);
      else
        // Break at false
        return false;
    }

    /// The extension list call implementation
    /// - it calls 'finalize()' on the current entry of the tuple
    /// - then broadcasts the extension call to the remaining tuple
    template <typename T, typename stepper_state_t>
    static bool
    finalize(const T& obs_tuple, stepper_state_t& state, const double h)
    {
      auto& this_extension = std::get<first>(obs_tuple);
      // Continue as long as evaluations are 'true'
      if (this_extension.finalize(state, h))
        return extension_list_impl<others...>::finalize(obs_tuple, state, h);
      else
        // Break at false
        return false;
    }
  };

  /// The extension list call implementation
  /// - it calls 'k()' on the last entry of the tuple
  template <typename last>
  struct extension_list_impl<last>
  {
    template <typename T, typename stepper_state_t>
    static bool
    k(T&                     obs_tuple,
      const stepper_state_t& state,
      Vector3D&              knew,
      const Vector3D&        bField,
      const int              i     = 0,
      const double           h     = 0,
      const Vector3D&        kprev = Vector3D())
    {
      auto& this_extension = std::get<last>(obs_tuple);
      return this_extension.k(state, knew, bField, i, h, kprev);
    }

    /// The extension list call implementation
    /// - it calls 'finalize()' on the last entry of the tuple
    template <typename T, typename stepper_state_t, typename stepper_data_t>
    static bool
    finalize(const T&              obs_tuple,
             stepper_state_t&      state,
             const double          h,
             const stepper_data_t& data,
             ActsMatrixD<7, 7>& D)
    {
      auto& this_extension = std::get<last>(obs_tuple);
      return this_extension.finalize(state, h, data, D);
    }

    /// The extension list call implementation
    /// - it calls 'finalize()' on the last entry of the tuple
    template <typename T, typename stepper_state_t>
    static bool
    finalize(const T& obs_tuple, stepper_state_t& state, const double h)
    {
      auto& this_extension = std::get<last>(obs_tuple);
      return this_extension.finalize(state, h);
    }
  };

  /// The empty extension list call implementation
  template <>
  struct extension_list_impl<>
  {
    template <typename T, typename stepper_state_t>
    static bool
    k(T& /*unused*/,
      const stepper_state_t& /*unused*/,
      Vector3D& /*unused*/,
      const Vector3D& /*unused*/,
      const int /*unused*/,
      const double /*unused*/,
      const Vector3D& /*unused*/)
    {
      return true;
    }

    /// The empty extension list call implementation
    template <typename T, typename stepper_state_t, typename stepper_data_t>
    static bool
    finalize(T& /*unused*/,
             stepper_state_t& /*unused*/,
             const double /*unused*/,
             const stepper_data_t& /*unused*/,
             ActsMatrixD<7, 7>& /*unused*/)
    {
      return true;
    }

    /// The empty extension list call implementation
    template <typename T, typename stepper_state_t>
    static bool
    finalize(T& /*unused*/,
             stepper_state_t& /*unused*/,
             const double /*unused*/)
    {
      return true;
    }
  };

}  // namespace detail
}  // namespace Acts

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
  namespace {

    struct extension_caller
    {
      template <typename extension, typename stepper_state_t>
      static void
      evaluatek1(extension& ext, stepper_state_t& state, const Vector3D& bField, Vector3D& k1)
      {
        //~ ext.k1(state, result.template get<detail::result_type_t<actor>>());
        ext.evaluatek1(state, bField, k1);
      }
    };
  }  // end of anonymous namespace

  /// The dummy list call implementation
  template <typename... extensions>
  struct extension_list_impl;

  /// The extension list call implementation
  /// - it calls 'extension' on the current entry of the tuple
  /// - then broadcasts the extension call to the remaining tuple
  template <typename first, typename... others>
  struct extension_list_impl<first, others...>
  {
    template <typename T, typename stepper_state_t>
    static void
    evaluatek1(T& obs_tuple, stepper_state_t& state, const Vector3D& bField, Vector3D& k1)
    {
      auto&    this_extension = std::get<first>(obs_tuple); // Extract extension
      this_extension.evaluatek1(state, bField, k1); // Call function
      extension_list_impl<others...>::evaluatek1(obs_tuple, state, bField, k1); // Forward to next element
    }
  };

  /// The extension list call implementation
  /// - it calls 'extension' on the last entry of the tuple
  template <typename last>
  struct extension_list_impl<last>
  {
    template <typename T, typename stepper_state_t>
    static void
    evaluatek1(T& obs_tuple, stepper_state_t& state, const Vector3D& bField, Vector3D& k1)
    {
      auto&    this_extension = std::get<last>(obs_tuple);
      extension_caller::evaluatek1(this_extension, state, bField, k1);
    }
  };

  /// The empty extension list call implementation
  template <>
  struct extension_list_impl<>
  {
    template <typename T, typename stepper_state_t>
    static void
    evaluatek1(T& /*unused*/,
           stepper_state_t& /*unused*/,
           const Vector3D& /*unused*/,
           Vector3D& /*unused*/)
    {
    }
  };

}  // namespace detail
}  // namespace Acts

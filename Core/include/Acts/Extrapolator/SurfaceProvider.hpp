// This file is part of the Acts project.
//
// Copyright (C) 2016-2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <sstream>
#include "Acts/Surfaces/Surface.hpp"

namespace Acts {

/// @brief User surface provider which can set the user
/// (or external surfaces) to the Navigator
///
/// This is implemented as an actor and should be provided
/// to the Propagator::Options<ActionList,AbortList> as one
/// of the abort list otions (ideally the first)
///
/// It will then set the user surfaces in the navigation state
/// and set itself to silent from then on.
///
/// @note This works only with the Acts::Navigator (not with)
/// the Acts::VoidNavigator.
struct SurfaceProvider
{

  /// Simple result struct to be returned
  /// It has all the SurfaceHit objects that
  /// are collected (and thus have been selected)
  struct this_result
  {
    bool initialized = false;
  };

  using result_type = this_result;

  /// Surface provider action for the ActionList of the Propagator
  ///
  /// @tparam propagator_state_t is the type of Propagator state
  ///
  /// @param[in,out] state is the mutable stepper state object
  /// @param[in,out] result is the mutable result object
  template <typename propagator_state_t>
  void
  operator()(propagator_state_t& state, result_type& result) const
  {
  }

  /// Pure observer interface
  /// - this does not apply to the surface collector
  template <typename propagator_state_t>
  void
  operator()(propagator_state_t& /*state*/) const
  {
  }
};

}  // namespace Acts

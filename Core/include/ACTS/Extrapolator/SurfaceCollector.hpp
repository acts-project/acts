// This file is part of the ACTS project.
//
// Copyright (C) 2016-2018 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <sstream>
#include "Acts/Surfaces/Surface.hpp"

namespace Acts {

/// The information to be writtern out per hit surface
struct SurfaceHit
{

  const Surface* surface = nullptr;
  Vector3D       position;
  Vector3D       direction;
};

/// A Surface Collector struct
/// templated with a Selector type
///
/// Whenever a surface is passed in the propagation
/// that satisfies the selector, it is recorded
/// for further usage in the flow.
template <typename Selector>
struct SurfaceCollector
{

  /// The selector used for this surface
  Selector this_selector;

  /// Simple result struct to be returned
  /// It has all the SurfaceHit objects that
  /// are collected (and thus have been selected)
  struct this_result
  {
    std::vector<SurfaceHit> collected;
  };

  typedef this_result result_type;

  /// Collector action for the ActionList of the Propagator
  /// It checks if the propagator state has a current surface,
  /// in which case the action is performed:
  /// - it records the surface given the configuration
  ///
  /// @tparam propagator_state_t is the type of Propagator state
  ///
  /// @param[in,out] state is the mutable stepper state object
  /// @param[in,out] result is the mutable result object
  template <typename propagator_state_t>
  void
  operator()(propagator_state_t& state, result_type& result) const
  {
    // a current surface has been assigned by the navigator
    //
    if (state.navigation.currentSurface
        && this_selector(*state.navigation.currentSurface)) {
      // create for recording
      SurfaceHit surface_hit;
      surface_hit.surface   = state.navigation.currentSurface;
      surface_hit.position  = state.stepping.position();
      surface_hit.direction = state.stepping.direction();
      // save if in the result
      result.collected.push_back(surface_hit);
    }
  }

  /// Pure observer interface
  /// - this does not apply to the surface collector
  template <typename propagator_state_t>
  void
  operator()(propagator_state_t&) const
  {
  }
};

}  // namespace Acts

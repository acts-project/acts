// This file is part of the ACTS project.
//
// Copyright (C) 2016-2018 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef ACTS_SURFACECOLLECTOR_H
#define ACTS_SURFACECOLLECTOR_H

#include <sstream>
#include "ACTS/Surfaces/Surface.hpp"

namespace Acts {

/// The information to written out per hit surface
struct SurfaceHit
{

  const Surface* surface = nullptr;
  Vector3D       position;
  Vector3D       direction;
};

/// A Surface Collector struct struct
/// templated with a Selector
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
  /// It checks if the cache has a current surface,
  /// in which case the action is performed:
  /// - it records the surface given the configuration
  ///
  /// @tparam propagator_cache_t is the type of Propagator cache
  /// @tparam stepper_cache_t is the type of Stepper cache
  ///
  /// @param pCache is the mutable stepper cache object
  /// @param sCache is the mutable stepper cache object
  /// @param result is the mutable result cache object
  template <typename propagator_cache_t, typename stepper_cache_t>
  void
  operator()(propagator_cache_t& pCache,
             stepper_cache_t&    sCache,
             result_type&        result) const
  {
    // a current surface has been assigned by the navigator
    //
    if (pCache.currentSurface && this_selector(*pCache.currentSurface)) {
      // create for recording
      SurfaceHit surface_hit;
      surface_hit.surface   = pCache.currentSurface;
      surface_hit.position  = sCache.position();
      surface_hit.direction = sCache.direction();
      // save if in the result
      result.collected.push_back(surface_hit);
    }
  }

  /// Pure observer interface
  /// - this does not apply to the surface collector
  template <typename propagator_cache_t, typename stepper_cache_t>
  void
  operator()(propagator_cache_t& pCache, stepper_cache_t& sCache) const
  {
    (void)pCache;
    (void)sCache;
  }
};
}

#endif

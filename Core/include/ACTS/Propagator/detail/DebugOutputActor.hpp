// This file is part of the ACTS project.
//
// Copyright (C) 2018 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef ACTS_DEBUG_OUTPUT_ACTOR_HPP
#define ACTS_DEBUG_OUTPUT_ACTOR_HPP

namespace Acts {

namespace detail {

  /// This is an actor that deals with the
  struct DebugOutputActor
  {

    /// mute the thing if you don't want any action
    bool mute = false;

    /// Simple result struct to be returned
    /// It collects the debug output string from the cache
    /// into which all actors and aborters can write
    struct this_result
    {
      std::string debugString = "";
    };

    typedef this_result result_type;

    /// Debug output action for the ActionList of the Propagator
    ///
    /// @tparam propagator_cache_t is the type of the Propagator cache
    /// @tparam stepper_cache_t is the type of Stepper cache,
    /// it is not used in this stepper
    ///
    /// @param pCache is the mutable propagator cache object
    /// @param result is the mutable result cache object
    template <typename propagator_cache_t, typename stepper_cache_t>
    void
    operator()(propagator_cache_t& pCache,
               stepper_cache_t&,
               result_type& result) const
    {
      // move the debug output from the cache to
      // to the output actor if it is not set to mute
      // only when the target is reached (or later otherwise triggered)
      if (!mute && (pCache.targetReached || pCache.navigationBreak))
        result.debugString = pCache.debugString;
    }

    /// Pure observer interface
    /// - this does not apply to the output collector
    template <typename propagator_cache_t, typename stepper_cache_t>
    void
    operator()(propagator_cache_t& pCache, stepper_cache_t& sCache) const
    {
      (void)pCache;
      (void)sCache;
    }
  };

}  // namespace detail
}  // namespace Acts

#endif

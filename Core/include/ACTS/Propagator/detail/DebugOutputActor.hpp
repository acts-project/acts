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
    struct this_result
    {
      std::string debugString = "";
    };

    typedef this_result result_type;

    /// Debug output action for the ActionList of the Propagator
    ///
    /// @tparam cache_t is the type of Stepper cache
    ///
    /// @param cache is the mutable stepper cache object
    /// @param result is the mutable result cache object
    template <typename cache_t>
    void
    operator()(cache_t& cache, result_type& result) const
    {
      // move the debug output from the cache to
      // to the output actor if it is not set to mute
      // only when the target is reached (or later otherwise triggered)
      if (!mute && cache.targetReached) result.debugString = cache.debugString;
    }

    /// Pure observer interface
    /// - this does not apply to the output collector
    template <typename cache_t>
    void
    operator()(cache_t& cache) const
    {
      (void)cache;
    }
  };

}  // namespace detail
}  // namespace Acts

#endif

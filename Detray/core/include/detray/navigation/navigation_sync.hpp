// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// Project include(s)
#include "detray/definitions/detail/qualifiers.hpp"

namespace detray::navigation {

/// @brief Policy that decides when the propagation loop runs the local
/// navigations that the navigator requested.
///
/// The propagation loop calls the policy once per iteration for every track,
/// whether or not the track requested a local navigation. The policy returns
/// whether the pending local navigations are run in this iteration. Until
/// then, a track with a pending request skips its actors and stepper.
///
/// This implementation never delays a local navigation.
struct no_sync {
  /// @param is_pending whether this track requested a local navigation
  /// @param n_waited number of iterations this track has been waiting
  DETRAY_HOST_DEVICE constexpr bool operator()(
      const bool /*is_pending*/, const unsigned int /*n_waited*/) const {
    return true;
  }
};

/// @brief Delays the local navigations until enough tracks in a warp request
/// one, so that the (very large) local navigation is run by many threads at
/// once.
///
/// The pending local navigations are run if at least @c min_pending threads
/// of the warp request one, if every active thread of the warp requests one,
/// or if one of the requesting threads has waited for @c max_wait iterations.
/// Outside of CUDA device code, the local navigations are never delayed.
struct warp_sync {
  /// Number of pending requests in the warp that triggers the navigations
  unsigned int min_pending{24u};
  /// Maximum number of iterations a track waits for its local navigation
  unsigned int max_wait{20u};

  /// @param is_pending whether this track requested a local navigation
  /// @param n_waited number of iterations this track has been waiting
  DETRAY_HOST_DEVICE inline bool operator()(
      [[maybe_unused]] const bool is_pending,
      [[maybe_unused]] const unsigned int n_waited) const {
#if defined(__CUDA_ARCH__)
    const unsigned int active = __activemask();
    const unsigned int pending = __ballot_sync(active, is_pending);
    const bool starving =
        __any_sync(active, is_pending && (n_waited >= max_wait));

    return (__popc(pending) >= static_cast<int>(min_pending)) ||
           (pending == active) || starving;
#else
    return true;
#endif
  }
};

}  // namespace detray::navigation

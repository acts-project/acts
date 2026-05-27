/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "traccc/definitions/qualifiers.hpp"
#include "traccc/utils/logging.hpp"

// Detray include(s).
#include <detray/propagator/base_actor.hpp>

// VecMem include(s).
#include <vecmem/utils/debug.hpp>

namespace traccc {

/// Aborter making sure that propagation would not exceed a certain step count
///
/// It is used mostly as a failsafe during the fitting stage. As long as
/// everything works well, this aborter should never trigger. It is here to
/// serve as a failsafe in case of a bug in the track finding / fitting
/// algorithm.
///
struct kalman_step_aborter : public detray::base_actor {

    /// The state of the aborter
    struct state {
        /// Maximum step count that a track can take to reach the next surface
        unsigned int max_steps = 100u;
        /// The current step count
        unsigned int step = 0u;
    };

    /// The operator to be called by the propagator
    ///
    /// @tparam propagator_state_t The type of the propagator state
    /// @param abrt_state The state of the aborter
    /// @param prop_state The state of the propagator
    ///
    template <typename propagator_state_t>
    TRACCC_HOST_DEVICE void operator()(state& abrt_state,
                                       propagator_state_t& prop_state) const {

        // Convenience reference to the navigation state.
        auto& navigation = prop_state.navigation();

        TRACCC_VERBOSE_HOST_DEVICE("Check Kalman step aborter...");

        // Reset the step count if the track is on a sensitive surface.
        if (navigation.is_on_sensitive()) {
            abrt_state.step = 0u;
        }

        // Abort if the step count exceeds the maximum allowed
        if (++(abrt_state.step) > abrt_state.max_steps) {
            navigation.abort(
                "Kalman Fitter: Maximum number of steps to reach next "
                "sensitive surface exceeded");
            prop_state.heartbeat(false);
        }
    }

};  // struct kalman_step_aborter

}  // namespace traccc

/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "traccc/definitions/primitives.hpp"
#include "traccc/definitions/qualifiers.hpp"
#include "traccc/utils/logging.hpp"

// detray include(s)
#include <detray/propagator/base_actor.hpp>

// System include(s)
#include <limits>

namespace traccc {

/// Aborter triggered when the next surface is reached
struct ckf_aborter : detray::base_actor {
    struct state {
        // minimal step length to prevent from staying on the same surface
        scalar min_step_length = 0.5f;
        /// Maximum step counts that track can make to reach the next surface
        unsigned int max_count = 100;

        bool success = false;
        unsigned int count = 0;

        scalar path_from_surface{0.f};
    };

    template <typename propagator_state_t>
    TRACCC_HOST_DEVICE void operator()(state &abrt_state,
                                       propagator_state_t &prop_state) const {

        auto &navigation = prop_state.navigation();
        const auto &stepping = prop_state.stepping();

        abrt_state.count++;
        abrt_state.path_from_surface += stepping.step_size();

        TRACCC_VERBOSE_HOST_DEVICE("Aborter: Checking CKF aborter");

        // Stop at the next sensitive surface
        if (navigation.is_on_sensitive() &&
            abrt_state.path_from_surface > abrt_state.min_step_length) {
            navigation.pause();
            prop_state.heartbeat(false);
            abrt_state.success = true;
            abrt_state.path_from_surface = 0.f;

            TRACCC_VERBOSE_HOST_DEVICE(
                "-> Found sensitive surface: %d",
                navigation.geometry_identifier().index());
        }

        if (abrt_state.count > abrt_state.max_count) {
            navigation.abort(
                "CKF: Maximum number of steps to reach next sensitive surface "
                "exceeded");
            prop_state.heartbeat(false);
        }
    }
};

}  // namespace traccc

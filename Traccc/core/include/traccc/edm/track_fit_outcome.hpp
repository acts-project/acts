/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2025 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// System include(s).
#include <cstdint>

namespace traccc {

/// Possible outcomes of a track fit
enum class track_fit_outcome : std::uint16_t {
    UNKNOWN,
    SUCCESS,
    FAILURE_NON_POSITIVE_NDF,
    FAILURE_NOT_ALL_FITTED,
    FAILURE_NOT_ALL_SMOOTHED,
    FAILURE_FITTER,
    FAILURE_SMOOTHER,
    FAILURE_FORWARD_PROPAGATION,
    FAILURE_BACKWARD_PROPAGATION,
    MAX_OUTCOME
};

}  // namespace traccc

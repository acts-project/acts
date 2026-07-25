/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2025 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Local include(s).
#include "traccc/fitting/kalman_filter/kalman_fitter.hpp"

// Detray include(s).
#include <detray/navigation/caching_navigator.hpp>
#include <detray/propagator/constrained_step.hpp>
#include <detray/propagator/rk_stepper.hpp>

// System include(s).
#include <type_traits>

namespace traccc::details {

/// Kalman fitter type for a specific detector and magnetic field type
///
/// @tparam detector_t The detector type to use
/// @tparam bfield_t   The magnetic field type to use
///
template <typename detector_t, typename bfield_t>
using kalman_fitter_t = kalman_fitter<
    detray::rk_stepper<
        bfield_t, typename detector_t::algebra_type,
        detray::constrained_step<traccc::scalar>,
        detray::stepper_rk_policy<traccc::scalar>,
        detray::stepping::void_inspector,
        static_cast<std::uint32_t>(
            detray::rk_stepper_flags::e_allow_covariance_transport)>,
    detray::caching_navigator<std::add_const_t<detector_t>>>;

}  // namespace traccc::details

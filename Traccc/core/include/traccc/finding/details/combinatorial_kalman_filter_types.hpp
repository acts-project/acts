/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2025 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Local include(s).
#include "traccc/definitions/common.hpp"
#include "traccc/definitions/primitives.hpp"
#include "traccc/finding/actors/ckf_aborter.hpp"

// Detray include(s).
#include "traccc/utils/propagation.hpp"

// System include(s).
#include <type_traits>

namespace traccc::details {

/// Stepper used in the Combinatorial Kalman Filter (CKF)
///
/// @tparam bfield_t The type of magnetic field to use
///
template <typename bfield_t>
using ckf_stepper_t = detray::rk_stepper<
    bfield_t, traccc::default_algebra, detray::constrained_step<traccc::scalar>,
    detray::stepper_rk_policy<traccc::scalar>, detray::stepping::void_inspector,
    static_cast<std::uint32_t>(
        detray::rk_stepper_flags::e_allow_covariance_transport)>;

/// Interactor used in the Combinatorial Kalman Filter (CKF)
using ckf_interactor_t =
    detray::actor::pointwise_material_interactor<traccc::default_algebra>;

/// Actor chain used in the Combinatorial Kalman Filter (CKF)
using ckf_actor_chain_t = detray::actor_chain<
    detray::actor::pathlimit_aborter<traccc::scalar>,
    detray::actor::parameter_updater<traccc::default_algebra, ckf_interactor_t>,
    detray::actor::momentum_aborter<traccc::scalar>, ckf_aborter>;

/// Propagator type used in the Combinatorial Kalman Filter (CKF)
///
/// @tparam detector_t The detector type to use
///
template <typename detector_t, typename bfield_t>
using ckf_propagator_t =
    detray::propagator<ckf_stepper_t<bfield_t>,
                       detray::caching_navigator<std::add_const_t<detector_t>>,
                       ckf_actor_chain_t>;

}  // namespace traccc::details

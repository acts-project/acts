/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Local include(s).
#include "traccc/definitions/common.hpp"
#include "traccc/definitions/primitives.hpp"
#include "traccc/finding/actors/measurement_kalman_updater.hpp"
#include "traccc/utils/propagation.hpp"

// System include(s).
#include <type_traits>

namespace traccc::details {

/// Stepper used in the Combinatorial Kalman Filter (CKF)
///
/// @tparam bfield_t The type of magnetic field to use
///
template <typename bfield_t>
using pkf_stepper_t = detray::rk_stepper<bfield_t, traccc::default_algebra>;

using parameter_updater_t = detray::actor::parameter_updater<
    traccc::default_algebra,
    detray::actor::pointwise_material_interactor<traccc::default_algebra>,
    traccc::measurement_updater<traccc::default_algebra>>;

template <typename surface_t>
using pkf_actor_chain_t =
    detray::actor_chain<detray::actor::pathlimit_aborter<traccc::scalar>,
                        detray::actor::surface_sequencer<surface_t>,
                        parameter_updater_t,
                        detray::actor::momentum_aborter<traccc::scalar>>;

template <typename detector_t, typename bfield_t>
using pkf_propagator_t =
    detray::propagator<pkf_stepper_t<bfield_t>,
                       detray::caching_navigator<std::add_const_t<detector_t>>,
                       pkf_actor_chain_t<typename detector_t::surface_type>>;

}  // namespace traccc::details

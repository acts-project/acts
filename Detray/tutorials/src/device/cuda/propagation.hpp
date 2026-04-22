// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// Project include(s).
#include "detray/core/detector.hpp"
#include "detray/definitions/units.hpp"
#include "detray/navigation/caching_navigator.hpp"
#include "detray/propagator/actors.hpp"
#include "detray/propagator/propagator.hpp"
#include "detray/propagator/rk_stepper.hpp"
#include "detray/tracks/tracks.hpp"

// Detray test include(s)
#include "detray/test/device/cuda/bfield.hpp"

// Tutorial include(s)
#include "detray/tutorial/types.hpp"

namespace detray::tutorial {

// Detector
using metadata_t = detray::tutorial::toy_metadata;
using detector_host_t = detector<metadata_t, host_container_types>;
using detector_device_t = detector<metadata_t, device_container_types>;

using algebra_t = metadata_t::algebra_type;
using scalar = detray::tutorial::scalar;

// Navigator
using navigator_t = caching_navigator<detector_device_t>;

// Stepper
using host_field_t = covfie::field<detray::bfield::inhom_bknd_t<scalar>>;
using device_field_t =
    covfie::field<detray::bfield::cuda::inhom_bknd_t<scalar>>;
using stepper_t = rk_stepper<device_field_t::view_t, algebra_t>;

// Actors
using actor_chain_t = actor_chain<
    actor::pathlimit_aborter<scalar>,
    actor::parameter_updater<algebra_t,
                             actor::pointwise_material_interactor<algebra_t>>>;

// Propagator
using propagator_t = propagator<stepper_t, navigator_t, actor_chain_t>;

// Free track parameters
using track_t = detray::free_track_parameters<detray::tutorial::algebra_t>;

/// Propagation tutorial function
void propagation(
    typename detector_host_t::view_type det_data,
    typename device_field_t::view_t field_data,
    const vecmem::data::vector_view<free_track_parameters<algebra_t>>
        tracks_data);

}  // namespace detray::tutorial

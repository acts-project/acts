/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2023-2025 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Local include(s).
#include "traccc/definitions/common.hpp"
#include "traccc/definitions/primitives.hpp"
#include "traccc/finding/measurement_selector.hpp"
#include "traccc/utils/particle.hpp"

// detray include(s).
#include <detray/propagator/propagation_config.hpp>

namespace traccc {

/// Configuration struct for track fitting
struct fitting_config {

    std::size_t n_iterations = 1;

    /// Propagation configuration
    detray::propagation::config propagation{};
    /// Measurement calibration configuration
    measurement_selector::config meas_calibration{};

    /// Minimum momentum for reconstructed tracks
    float min_p = 100.f * traccc::unit<float>::MeV;
    float min_pT = 600.f * traccc::unit<float>::MeV;

    /// Particle hypothesis
    traccc::pdg_particle<traccc::scalar> ptc_hypothesis =
        traccc::pion_plus<traccc::scalar>();

    /// Smoothing with backward filter
    traccc::scalar covariance_inflation_factor = 1e3f;
    unsigned int surface_sequence_size_factor = 5u;
    unsigned int min_surface_sequence_capacity = 100u;
    bool do_precise_hole_count = false;
};

}  // namespace traccc

/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// traccc include(s).
#include "traccc/definitions/primitives.hpp"
#include "traccc/definitions/qualifiers.hpp"
#include "traccc/edm/track_parameters.hpp"

// detray include(s).
#include <detray/definitions/pdg_particle.hpp>

// System include(s).
#include <stdexcept>

namespace traccc {

template <typename scalar_t>
using pdg_particle = detray::pdg_particle<scalar_t>;

template <typename scalar_t>
using electron = detray::electron<scalar_t>;
template <typename scalar_t>
using positron = detray::positron<scalar_t>;
template <typename scalar_t>
using muon = detray::muon<scalar_t>;
template <typename scalar_t>
using antimuon = detray::antimuon<scalar_t>;
template <typename scalar_t>
using pion_plus = detray::pion_plus<scalar_t>;
template <typename scalar_t>
using pion_minus = detray::pion_minus<scalar_t>;

namespace detail {

template <typename scalar_t>
TRACCC_HOST_DEVICE inline traccc::pdg_particle<scalar_t>
particle_from_pdg_number(const int pdg_num) {

    switch (pdg_num) {
        case 11:
            return detray::electron<scalar_t>();
        case -11:
            return detray::positron<scalar_t>();
        case 13:
            return detray::muon<scalar_t>();
        case -13:
            return detray::antimuon<scalar_t>();
        case 211:
            return detray::pion_plus<scalar_t>();
        case -211:
            return detray::pion_minus<scalar_t>();
    }

    // TODO: Replace with `detray::invalid` in the future
    return traccc::pdg_particle<scalar_t>(0, 0.f, 0.f);
}

// Apply the charge operator to return the antimatter
template <typename scalar_t>
TRACCC_HOST_DEVICE inline traccc::pdg_particle<scalar_t> charge_conjugation(
    const traccc::pdg_particle<scalar_t>& ptc) {

    const auto pdg_num = ptc.pdg_num();

    switch (pdg_num) {
        case 11:
            return detray::positron<scalar_t>();
        case -11:
            return detray::electron<scalar_t>();
        case 13:
            return detray::antimuon<scalar_t>();
        case -13:
            return detray::muon<scalar_t>();
        case 211:
            return detray::pion_minus<scalar_t>();
        case -211:
            return detray::pion_plus<scalar_t>();
    }

    // TODO: Replace with `detray::invalid` in the future
    return traccc::pdg_particle<scalar_t>(0, 0.f, 0.f);
}

// Return the consistent particle type based on the particle hypothesis and the
// charge of the track parameters
template <typename scalar_t>
TRACCC_HOST_DEVICE inline traccc::pdg_particle<scalar_t>
correct_particle_hypothesis(
    const traccc::pdg_particle<scalar_t>& ptc_hypothesis,
    const bound_track_parameters<>& params) {

    if (ptc_hypothesis.charge() * params.qop() > 0.f) {
        return ptc_hypothesis;
    } else {
        return charge_conjugation(ptc_hypothesis);
    }
}

template <typename scalar_t>
TRACCC_HOST_DEVICE inline traccc::pdg_particle<scalar_t>
correct_particle_hypothesis(
    const traccc::pdg_particle<scalar_t>& ptc_hypothesis,
    const free_track_parameters<>& params) {

    if (ptc_hypothesis.charge() * params.qop() > 0.f) {
        return ptc_hypothesis;
    } else {
        return charge_conjugation(ptc_hypothesis);
    }
}

}  // namespace detail

}  // namespace traccc

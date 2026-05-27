/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Library include(s).
#include "traccc/edm/track_parameters.hpp"
#include "traccc/utils/logging.hpp"
#include "traccc/utils/particle.hpp"

// detray include(s).
#include <detray/definitions/navigation.hpp>  // < navigation::direction
#include <detray/geometry/identifier.hpp>
#include <detray/geometry/tracking_surface.hpp>
#include <detray/propagator/actors/pointwise_material_interactor.hpp>

// System include(s).
#include <random>

namespace traccc {

/// Seed track parameter generator
template <typename detector_t>
struct seed_generator {
    using algebra_type = typename detector_t::algebra_type;
    using ctx_t = typename detector_t::geometry_context;

    /// Constructor with detector
    ///
    /// @param det input detector
    /// @param stddevs standard deviations for parameter smearing
    seed_generator(const detector_t& det,
                   const std::array<scalar, e_bound_size>& stddevs,
                   const std::size_t sd = 0, ctx_t ctx = {})
        : m_detector(det), m_stddevs(stddevs), m_ctx(ctx) {
        m_generator.seed(static_cast<std::mt19937::result_type>(sd));
    }

    /// Seed generator operation
    ///
    /// @param vertex vertex of particle
    /// @param stddevs standard deviations for track parameter smearing
    bound_track_parameters<algebra_type> operator()(
        const detray::geometry::identifier surface_link,
        const free_track_parameters<algebra_type>& free_param,
        const traccc::pdg_particle<scalar>& ptc_type) {

        // Get bound parameter
        const detray::tracking_surface sf{m_detector, surface_link};

        TRACCC_DEBUG_HOST("Surface:\n" << sf);
        TRACCC_DEBUG_HOST("Input free param. on surface:\n" << free_param);

        auto bound_vec = sf.free_to_bound_vector(m_ctx, free_param);
        auto bound_cov = matrix::identity<traccc::bound_matrix<algebra_type>>();

        bound_track_parameters<algebra_type> bound_param{surface_link,
                                                         bound_vec, bound_cov};

        TRACCC_DEBUG_HOST("-> Bound param.:\n" << bound_param);

        // Type definitions
        using interactor_type =
            detray::actor::pointwise_material_interactor<algebra_type>;

        assert(ptc_type.charge() * bound_param.qop() > 0.f);

        // Apply interactor
        TRACCC_DEBUG_HOST("Update material:");
        if (sf.has_material()) {
            typename interactor_type::state interactor_state;
            interactor_state.do_multiple_scattering = false;
            interactor_type{}.update(
                m_ctx, ptc_type, bound_param, interactor_state,
                static_cast<int>(detray::navigation::direction::e_backward),
                sf);
        }

        for (std::size_t i = 0; i < e_bound_size; i++) {

            if (m_stddevs[i] != scalar{0}) {
                bound_param[i] = std::normal_distribution<scalar>(
                    bound_param[i], m_stddevs[i])(m_generator);
            }

            getter::element(bound_param.covariance(), i, i) =
                m_stddevs[i] * m_stddevs[i];
        }

        TRACCC_VERBOSE_HOST("Generated seed param.:\n" << bound_param);

        assert(!bound_param.is_invalid());

        return bound_param;
    }

    private:
    // Random generator
    std::random_device m_rd{};
    std::mt19937 m_generator{m_rd()};

    // Detector object
    const detector_t& m_detector;
    /// Standard deviations for parameter smearing
    std::array<scalar, e_bound_size> m_stddevs;
    ctx_t m_ctx;
};

}  // namespace traccc

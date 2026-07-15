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
#include "traccc/utils/trigonometric_helpers.hpp"

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

    /// Configure the seed generator
    struct config {
        // Smearing parameters
        /// Constant term of the loc0 resolution.
        double sigma_loc0 = 20 * traccc::unit<scalar>::um;
        /// Pt-dependent loc0 resolution of the form sigma_loc0 =
        /// A*exp(-1.*abs(B)*pt).
        double sigma_loc0_pT_a = 30 * traccc::unit<scalar>::um;
        double sigma_loc0_pT_b = 0.3 / traccc::unit<scalar>::GeV;
        /// Constant term of the loc1 resolution.
        double sigma_loc1 = 20 * traccc::unit<scalar>::um;
        /// Pt-dependent loc1 resolution of the form sigma_loc1 =
        /// A*exp(-1.*abs(B)*pt).
        double sigma_loc1_pT_a = 30 * traccc::unit<scalar>::um;
        double sigma_loc1_pT_b = 0.3 / traccc::unit<scalar>::GeV;
        /// Time resolution.
        double sigma_time = 1 * traccc::unit<scalar>::ns;
        /// Phi angular resolution.
        double sigma_phi = 1 * traccc::unit<scalar>::degree;
        /// Theta angular resolution.
        double sigma_theta = 1 * traccc::unit<scalar>::degree;
        /// Relative transverse momentum resolution.
        double sigma_pT_rel = 0.05;

        /// Optional. Initial sigmas for the track parameters which overwrites
        /// the smearing params if set.
        std::optional<std::array<double, e_bound_size>> initial_sigmas;
        /// Initial sigma(q/pt) for the track parameters.
        /// @note The resulting q/p sigma is added to the one in `initialSigmas`
        double initialsigma_qopt =
            0.1 * traccc::unit<scalar>::e / traccc::unit<scalar>::GeV;
        /// Initial sigma(pt)/pt for the track parameters.
        /// @note The resulting q/p sigma is added to the one in `initialSigmas`
        double initialsigma_pT_rel = 0.1;

        /// Inflation factors for the variances
        std::array<double, e_bound_size> cov_inflation{1., 1., 1., 1., 1., 1.};
    };

    /// Constructor with detector
    ///
    /// @param det input detector
    /// @param stddevs standard deviations for parameter smearing
    explicit seed_generator(const detector_t& det, const config& cfg = {},
                            const std::size_t sd = std::mt19937::default_seed,
                            ctx_t ctx = {})
        : m_detector(det), m_ctx(ctx), m_cfg(cfg) {
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

        // Call the smearing operator
        (*this)(bound_param, ptc_type);

        bound_param.set_surface_link(surface_link);

        return bound_param;
    }

    /// Seed generator operation
    ///
    /// @param vertex vertex of particle
    /// @param stddevs standard deviations for track parameter smearing
    void operator()(bound_track_parameters<algebra_type>& bound_param,
                    const traccc::pdg_particle<scalar>& ptc_type) {

        const scalar q{ptc_type.charge()};
        const auto pos = bound_param.bound_local();
        const auto time = bound_param.time();
        const auto phi = bound_param.phi();
        const auto theta = bound_param.theta();
        const auto pt = bound_param.pT(q);
        const auto qop = bound_param.qop();

        // Build the track covariance matrix using the smearing sigmas
        const auto sigmas = generate_smearing_sigmas(pt, qop, theta);

        // Smear the position/time
        // Note that we smear d0 and z0 in the perigee frame
        const double smeared_loc0 = std::normal_distribution<double>(
            pos[0], sigmas[e_bound_loc0])(m_generator);
        const double smeared_loc1 = std::normal_distribution<double>(
            pos[1], sigmas[e_bound_loc1])(m_generator);
        bound_param.set_bound_local({static_cast<scalar>(smeared_loc0),
                                     static_cast<scalar>(smeared_loc1)});

        // Time
        bound_param.set_time(
            static_cast<scalar>(std::normal_distribution<double>(
                time, sigmas[e_bound_time])(m_generator)));

        // Smear direction angles phi,theta ensuring correct bounds
        const double smeared_phi = std::normal_distribution<double>(
            phi, sigmas[e_bound_phi])(m_generator);
        const double smeared_theta = std::normal_distribution<double>(
            theta, sigmas[e_bound_theta])(m_generator);
        const auto [new_phi, new_theta] =
            detail::wrap_phi_theta(static_cast<scalar>(smeared_phi),
                                   static_cast<scalar>(smeared_theta));
        bound_param.set_phi(static_cast<scalar>(new_phi));
        bound_param.set_theta(static_cast<scalar>(new_theta));

        // Compute smeared q/p
        bound_param.set_qop(
            static_cast<scalar>(std::normal_distribution<double>(
                qop, sigmas[e_bound_qoverp])(m_generator)));

        // Generate the covariances
        generate_initial_covariance(bound_param, ptc_type);

        // We are not on a surface, but in the curvilinear frame
        bound_param.set_surface_link(detray::geometry::identifier{});
    }

    /// Generate an array of standard deviations
    ///
    /// @param pt transverse momentum of the particle
    /// @param qop q/p of the particle
    /// @param theta global theta angle of the particle direction
    ///
    /// @returns an array of stddeviations corresponding to each track parameter
    std::array<double, e_bound_size> generate_smearing_sigmas(
        const double pt, const double qop, const double theta) {

        if (m_cfg.initial_sigmas.has_value()) {
            return *m_cfg.initial_sigmas;
        } else {
            // Momentum dependent standard deviations
            const double sigma_loc_0 =
                m_cfg.sigma_loc0 +
                m_cfg.sigma_loc0_pT_a *
                    math::exp(-1.0 * math::fabs(m_cfg.sigma_loc0_pT_b) * pt);

            const double sigma_loc_1 =
                m_cfg.sigma_loc1 +
                m_cfg.sigma_loc1_pT_a *
                    math::exp(-1.0 * math::fabs(m_cfg.sigma_loc1_pT_b) * pt);

            const double sigma_qop = math::sqrt(
                math::pow(m_cfg.sigma_pT_rel * qop, 2) +
                math::pow(m_cfg.sigma_theta * (qop * math::tan(theta)), 2));

            return {sigma_loc_0,       sigma_loc_1, m_cfg.sigma_phi,
                    m_cfg.sigma_theta, sigma_qop,   m_cfg.sigma_time};
        }
    }

    /// Generate the covariance matrix for the track
    ///
    /// @param[in, out] bound_param the bound track parameters
    /// @param ptc_type the particle hypothesis
    void generate_initial_covariance(
        bound_track_parameters<algebra_type>& bound_param,
        const traccc::pdg_particle<scalar>& ptc_type) {

        const scalar q{ptc_type.charge()};

        // Build the track covariance matrix using the smearing sigmas
        const auto sigmas = generate_smearing_sigmas(
            bound_param.pT(q), bound_param.qop(), bound_param.theta());

        for (std::size_t i = e_bound_loc0; i < e_bound_size; ++i) {
            double sigma = sigmas[i];
            double variance = sigma * sigma;

            // Inflate the initial covariance
            variance *= m_cfg.cov_inflation[i];

            getter::element(bound_param.covariance(), i, i) =
                static_cast<scalar>(variance);
        }

        assert(!bound_param.is_invalid());
    }

    private:
    /// Random generator
    std::random_device m_rd{};
    std::mt19937 m_generator{m_rd()};

    /// Detector object
    const detector_t& m_detector;
    /// Geometry context
    ctx_t m_ctx;

    /// Seed generator configuration
    config m_cfg{};
};

}  // namespace traccc

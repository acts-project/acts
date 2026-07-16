/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "traccc/definitions/math.hpp"
#include "traccc/definitions/primitives.hpp"
#include "traccc/definitions/qualifiers.hpp"
#include "traccc/edm/measurement_collection.hpp"
#include "traccc/edm/measurement_helpers.hpp"
#include "traccc/edm/track_container.hpp"
#include "traccc/edm/track_parameters.hpp"
#include "traccc/edm/track_state_collection.hpp"
#include "traccc/fitting/fitting_config.hpp"
#include "traccc/utils/prob.hpp"

// detray include(s).
#include <detray/geometry/tracking_surface.hpp>
#include <detray/tracks/bound_track_parameters.hpp>
#include <detray/utils/geometry_utils.hpp>

// vecmem include(s)
#include <vecmem/containers/device_vector.hpp>

// System include(s).
#include <iostream>
#include <limits>
#include <string>

namespace traccc {

/// Triplet fitting algorithm to fit a single track

template <typename detector_t, typename bfield_t>
class triplet_fitter {

    public:
    // Algebra type
    using algebra_type = typename detector_t::algebra_type;

    // Configuration type
    using config_type = fitting_config;

    // Matrix types
    using size_type = detray::dindex_type<algebra_type>;
    template <size_type ROWS, size_type COLS>
    using matrix_type = detray::dmatrix<algebra_type, ROWS, COLS>;

    /// Constructor with a detector
    ///
    /// @param det the detector object
    /// @param field magnetic field
    /// @param cfg fitter configuration
    TRACCC_HOST_DEVICE
    triplet_fitter(const detector_t& det, const bfield_t& field,
                   const config_type& cfg)
        : m_detector(det), m_field(field), m_cfg(cfg) {}

    // Convenience function: wrap angle between -Pi and Pi
    static inline scalar wrap_pi_mpi(scalar angle) {
        // Map angle to [ -2*PI, 2*PI ]
        const scalar two_pi = 2.0f * static_cast<scalar>(M_PI);
        angle = std::fmod(angle + static_cast<scalar>(M_PI), two_pi);

        // Handle negative results from fmod to ensure we are in [ 0, 2*PI ]
        if (angle < 0) {
            angle += two_pi;
        }

        // Shift back to [ -PI, PI ]
        return angle - static_cast<scalar>(M_PI);
    }

    // Convenience function: difference of two azimuthal angles
    static inline scalar delta_phi(scalar phi_a, scalar phi_b) {
        scalar dphi =
            math::fmod(phi_a - phi_b, 2.f * static_cast<scalar>(M_PI));

        if (dphi > static_cast<scalar>(M_PI))
            dphi -= 2.f * static_cast<scalar>(M_PI);
        if (dphi < -1.f * static_cast<scalar>(M_PI))
            dphi += 2.f * static_cast<scalar>(M_PI);
        return dphi;
    }

    // Segment parameters
    struct segment_par {
        scalar phi;
        scalar theta;
        scalar rho_phi;
        scalar rho_theta;
        scalar nu;
        scalar tau;
    };

    // Hit gradients of segment parameters
    struct segment_grad {

        // Construct in-place for containers
        segment_grad(scalar grad_phi, scalar grad_theta, scalar grad_rho_phi,
                     scalar grad_rho_theta, scalar grad_nu, scalar grad_tau)
            : f_phi{grad_phi},
              f_theta{grad_theta},
              f_rho_phi{grad_rho_phi},
              f_rho_theta{grad_rho_theta},
              f_nu{grad_nu},
              f_tau{grad_tau} {}

        scalar f_phi;
        scalar f_theta;
        scalar f_rho_phi;
        scalar f_rho_theta;
        scalar f_nu;
        scalar f_tau;
    };

    // Triplet struct
    struct triplet {

        /// Default construct
        triplet() = default;

        /// Construct with three hits
        ///
        /// @param hit_0 position of first hit
        /// @param hit_1 position of second hit
        /// @param hit_2 position of third hit
        triplet(const point3& hit_0, const point3& hit_1, const point3& hit_2)
            : m_hit_pos{hit_0, hit_1, hit_2} {}

        /// Custom copy
        ///
        /// (only hit positions and measurements copied)
        triplet(const triplet& t)
            : m_hit_pos{t.m_hit_pos}, m_meas_idx{t.m_meas_idx} {}

        /// Shift a hit by measurement uncertainty
        ///
        /// @param hit_idx hit index
        /// @param dir_idx direction index
        /// @param multipler multipler
        /// @param detector detector
        void shiftHit(const unsigned& hit_idx, const unsigned& dir_idx,
                      const scalar& multipler,
                      edm::measurement_collection::const_device measurements,
                      const detector_t& detector) {
            // The shifts are applied to the global
            // positions of the hits

            // Get local positions and variances of measurement
            point2 loc_pos{};
            point2 loc_var{};
            const auto meas = measurements.at(m_meas_idx[hit_idx]);
            for (std::size_t dim : meas.subspace()) {
                loc_pos[dim] =
                    measurements.at(m_meas_idx[hit_idx]).local_position()[dim];
                loc_var[dim] =
                    measurements.at(m_meas_idx[hit_idx]).local_variance()[dim];
            }

            // Shift local position
            assert(dir_idx <=
                   measurements.at(m_meas_idx[hit_idx]).dimensions());
            loc_pos[dir_idx] += multipler * math::sqrt(loc_var[dir_idx]);

            // Surface
            detray::tracking_surface sf{
                detector, measurements.at(m_meas_idx[hit_idx]).surface_link()};

            // Convert shifted position to global coordinates
            point3 glob_pos = sf.local_to_global({}, loc_pos, {});

            // Save shifted position in triplet object
            this->m_hit_pos[hit_idx] = glob_pos;
        }

        // Hit positions
        std::array<point3, 3u> m_hit_pos;

        // Triplet parameters
        scalar m_phi_0{};
        scalar m_theta_0{};
        scalar m_rho_phi{};
        scalar m_rho_theta{};

        // Segment parameters
        segment_par m_seg_1;
        segment_par m_seg_2;

        scalar m_sigma_MS{};  // MS uncertainty
        scalar m_theta{};     // Polar angle

        // Triplet hit position derivatives
        // (h * sqrt(Dinv))
        // (per hit per uncertainty direction)
        vecmem::vector<scalar> m_hd_theta;
        vecmem::vector<scalar> m_hd_phi;

        // Triplet segment gradients
        vecmem::vector<segment_grad> m_fseg_1;
        vecmem::vector<segment_grad> m_fseg_2;

        // Measurements for getting position uncertainties (hit shifts)
        // and surface orientation (scattering estimation)
        // std::array<typename
        // edm::measurement_collection::host::proxy_type, 3u>
        // m_meas;

        // No default construction for measurement proxies
        // can I have an array of indices instead ?
        std::array<unsigned int, 3u> m_meas_idx;
    };

    /// Helper function - Make triplets
    ///
    /// Makes triplets from consecutive measurements on the track
    ///
    TRACCC_HOST_DEVICE
    void make_triplets(const vecmem::vector<unsigned int>& in_measurements,
                       edm::measurement_collection::const_device measurements) {

        // Assuming no holes
        const size_t n_triplets = in_measurements.size() - 2;

        // Clear triplets from last candidate
        m_triplets.clear();
        m_meas_sum_dims = 0;

        m_triplets.reserve(n_triplets);

        // loop over measurements (track states) in candidate
        for (size_t i = 0; i < n_triplets; ++i) {

            // Get track states (and measurements)
            auto meas_0 = measurements.at(in_measurements[i]);
            auto meas_1 = measurements.at(in_measurements[i + 1]);
            auto meas_2 = measurements.at(in_measurements[i + 2]);

            // Get surfaces
            detray::tracking_surface meas_0_sf(m_detector,
                                               meas_0.surface_link());
            detray::tracking_surface meas_1_sf(m_detector,
                                               meas_1.surface_link());
            detray::tracking_surface meas_2_sf(m_detector,
                                               meas_2.surface_link());

            // Get measurement local positions
            point2 loc_2d_0{};
            for (std::size_t dim : meas_0.subspace()) {
                loc_2d_0[dim] = meas_0.local_position()[dim];
            }

            point2 loc_2d_1{};
            for (std::size_t dim : meas_1.subspace()) {
                loc_2d_1[dim] = meas_1.local_position()[dim];
            }

            point2 loc_2d_2{};
            for (std::size_t dim : meas_2.subspace()) {
                loc_2d_2[dim] = meas_2.local_position()[dim];
            }

            // Convert to global
            point3 glob_3d_0 = meas_0_sf.local_to_global({}, loc_2d_0, {});
            point3 glob_3d_1 = meas_1_sf.local_to_global({}, loc_2d_1, {});
            point3 glob_3d_2 = meas_2_sf.local_to_global({}, loc_2d_2, {});

            // Make triplet
            triplet t(glob_3d_0, glob_3d_1, glob_3d_2);

            t.m_meas_idx[0] = in_measurements[i];
            t.m_meas_idx[1] = in_measurements[i + 1];
            t.m_meas_idx[2] = in_measurements[i + 2];

            // keep track of dimensions
            m_meas_sum_dims += static_cast<scalar>(meas_0.dimensions());
            if (i == n_triplets - 1u) {
                m_meas_sum_dims += static_cast<scalar>(meas_1.dimensions());
                m_meas_sum_dims += static_cast<scalar>(meas_2.dimensions());
            }

            m_triplets.push_back(t);
        }
    }

    /// Estimate track direction at the scattering plane
    /// using hits 0, 2 and the circle solution
    ///
    /// @param t triplet
    /// @param c_perp curvature of circle
    /// @param dir track direction (output)
    void estimateTrackDirection(const triplet& t, const scalar& c_perp,
                                vector3& dir) {

        // Mid-point between hits 0 & 2
        vector2 m{0.5f * (t.m_hit_pos[0][0] + t.m_hit_pos[2][0]),
                  0.5f * (t.m_hit_pos[0][1] + t.m_hit_pos[2][1])};

        vector3 x_02{t.m_hit_pos[2] - t.m_hit_pos[0]};
        scalar d_02 = detray::algebra::array::perp(x_02);

        // Direction perpendicular to vector joining hits 0 & 2
        vector2 n{(t.m_hit_pos[2][1] - t.m_hit_pos[0][1]) / d_02,
                  (t.m_hit_pos[0][0] - t.m_hit_pos[2][0]) / d_02};

        scalar perp_d =
            math::sqrt(1.f / (c_perp * c_perp) - 0.25f * (d_02 * d_02));

        // Centre can be on either side of line joining hits 0 & 2
        std::array<vector2, 2u> c;

        c[0] = vector2{m[0] + n[0] * perp_d, m[1] + n[1] * perp_d};
        c[1] = vector2{m[0] - n[0] * perp_d, m[1] - n[1] * perp_d};

        vector2 x1{t.m_hit_pos[1][0], t.m_hit_pos[1][1]};

        // Choose the correct centre
        vector2 c_correct{0.f, 0.f};
        for (const vector2& c_i : c) {
            // Centre of the circle cannot be
            // on the same side of the line
            // connecting hits 0 & 2 as hit 1
            if (vector::dot(x1 - m, c_i - m) < 0.f) {
                c_correct = c_i;
                break;
            }
        }

        if (detray::algebra::array::norm(c_correct) == 0.f or
            detray::algebra::array::norm(x1 - m) == 0.f) {
            // Use vector joining hits
            // 0 and 2 as track direction if
            // center calculation fails or three
            // hits lie on a straight line
            dir = vector::normalize(x_02);
        }

        else {
            // Use circle solution to get track direction

            vector2 r1 = x1 - c_correct;
            vector2 tangent2D{r1[1], -1.f * r1[0]};

            // tangent direction along trajectory
            vector3 x_12{t.m_hit_pos[2] - t.m_hit_pos[1]};
            if (vector::dot(tangent2D, vector2{x_12[0], x_12[1]}) < 0.f)
                tangent2D = -1.f * tangent2D;

            vector2 tangent2D_norm = math::sin(t.m_theta) /
                                     detray::algebra::array::norm(tangent2D) *
                                     tangent2D;

            // track tangent normalized to 1
            dir[0] = tangent2D_norm[0];
            dir[1] = tangent2D_norm[1];
            dir[2] = math::cos(t.m_theta);
        }
    }

    /// Helper function - Linearize triplet
    ///
    /// Calculates triplet parameters by linearizing around circle solution
    /// Estimates triplet polar angle and multiple scattering uncertainty
    ///
    /// @param t Triplet to linearize
    ///
    TRACCC_HOST_DEVICE void linearize_triplet(
        triplet& t, edm::measurement_collection::const_device measurements) {

        // Vectors joining hits
        vector3 x_01{t.m_hit_pos[1] - t.m_hit_pos[0]};
        vector3 x_12{t.m_hit_pos[2] - t.m_hit_pos[1]};
        vector3 x_02{t.m_hit_pos[2] - t.m_hit_pos[0]};

        // Transverse distances
        scalar d_01 = detray::algebra::array::perp(x_01);
        scalar d_12 = detray::algebra::array::perp(x_12);
        scalar d_02 = detray::algebra::array::perp(x_02);

        // Longitudinal distances
        scalar z_01 = x_01[2];
        scalar z_12 = x_12[2];

        // Azimuthal angles of transverse vectors
        scalar phi_01 = math::atan2(x_01[1], x_01[0]);
        scalar phi_12 = math::atan2(x_12[1], x_12[0]);

        // Calculation of circle curvature and hence the entire
        // linearization will fail for very low (or 0) transverse
        // distances between hits. The default initialized (0)
        // values of triplet parameters are returned in this case.
        constexpr scalar d_transverse_lim = 10e-6f;

        // Curvature of circle in transverse plane
        scalar c_perp;
        if ((d_01 > d_transverse_lim and d_12 > d_transverse_lim and
             d_02 > d_transverse_lim)) {
            // TODO: x-prod evaluates to -ve, might have to be reversed
            c_perp = 2.f * math::fabs((vector::cross(x_01, x_12))[2]) /
                     (d_01 * d_12 * d_02);
        } else {
            return;
        }

        // Parameters of the arc segments //

        // Transverse bending angles
        // (with the assumption that smaller of the two
        // transverse bending angles has |PHI_xC| < PI)
        // 1, 2: segment
        // C: circle solution in transverse plane

        scalar sin_phiHalf_1C = 0.5f * d_01 * c_perp;
        scalar sin_phiHalf_2C = 0.5f * d_12 * c_perp;

        scalar phi_1C, phi_2C;

        phi_1C = 2.f * math::asin(sin_phiHalf_1C);
        phi_2C = 2.f * math::asin(sin_phiHalf_2C);

        // Bending parameters
        scalar Xphi_1C = 1.f;
        scalar Xphi_2C = 1.f;
        if (d_01 > 0.f)
            Xphi_1C = 0.5f * phi_1C / sin_phiHalf_1C;
        if (d_12 > 0.f)
            Xphi_2C = 0.5f * phi_2C / sin_phiHalf_2C;

        scalar cos_phiHalf_1C =
            math::sqrt(1.f - sin_phiHalf_1C * sin_phiHalf_1C);
        scalar cos_phiHalf_2C =
            math::sqrt(1.f - sin_phiHalf_2C * sin_phiHalf_2C);

        // Transverse arc lengths
        scalar t_1C, t_2C;

        if (c_perp == 0.f) {
            t_1C = d_01;
            t_2C = d_12;
        } else {
            t_1C = phi_1C / c_perp;
            t_2C = phi_2C / c_perp;
        }

        // Total arc lengths
        scalar s_1C_2 = z_01 * z_01 + t_1C * t_1C;
        scalar s_2C_2 = z_12 * z_12 + t_2C * t_2C;

        // Polar angles
        scalar theta_1C = math::atan2(t_1C, z_01);
        scalar theta_2C = math::atan2(t_2C, z_12);

        scalar cot_th_1C = z_01 / t_1C;
        scalar cot_th_2C = z_12 / t_2C;

        scalar cos_th_1C_2 = z_01 * z_01 / s_1C_2;
        scalar cos_th_2C_2 = z_12 * z_12 / s_2C_2;

        scalar sin_th_1C_2 = 1.f - cos_th_1C_2;
        scalar sin_th_2C_2 = 1.f - cos_th_2C_2;

        scalar sin_th_1C = math::sqrt(sin_th_1C_2);
        scalar sin_th_2C = math::sqrt(sin_th_2C_2);

        // Estimate polar angle of the triplet
        t.m_theta = 0.5f * (theta_1C + theta_2C);

        // Track direction at scattering surface
        vector3 tangent3D;
        estimateTrackDirection(t, c_perp, tangent3D);

        // Estimate MS-uncertainty
        // (track direction used here
        // to get precise thickness of material)

        detray::tracking_surface scat_sf(
            m_detector, measurements.at(t.m_meas_idx[1]).surface_link());

        // effective thickness
        scalar t_eff =
            mat_scatter /
            detray::cos_angle({}, scat_sf, tangent3D,
                              edm::get_measurement_local<algebra_type>(
                                  measurements.at(t.m_meas_idx[1])));

        auto scattering_unc = [](scalar curvature_3D, scalar eff_thickness,
                                 vector3 field_strength_vector) {
            return math::fabs(curvature_3D) * 45.0311528518f *
                   math::sqrt(eff_thickness) * unit<scalar>::T /
                   field_strength_vector[2] *
                   (1.f + 0.038f * math::log(eff_thickness));
        };

        const auto B_field =
            m_field.at(t.m_hit_pos[1][0], t.m_hit_pos[1][1], t.m_hit_pos[1][2]);
        vector3 B_vec;
        B_vec[0u] = B_field[0u];
        B_vec[1u] = B_field[1u];
        B_vec[2u] = B_field[2u];

        // 3D curvatures of segments
        scalar c_3D_1C = c_perp * sin_th_1C;
        scalar c_3D_2C = c_perp * sin_th_2C;

        scalar c3D_lin = 0.5f * (c_3D_1C + c_3D_2C);
        t.m_sigma_MS = scattering_unc(c3D_lin, t_eff, B_vec);

        // Index parameters
        scalar n_1C =
            1.f / (Xphi_1C * cos_phiHalf_1C * sin_th_1C_2 + cos_th_1C_2);
        scalar n_2C =
            1.f / (Xphi_2C * cos_phiHalf_2C * sin_th_2C_2 + cos_th_2C_2);

        // Segment parameters

        // rho_phi
        t.m_seg_1.rho_phi = t_1C * n_1C / sin_th_1C;
        t.m_seg_2.rho_phi = t_2C * n_2C / sin_th_2C;

        // tau
        t.m_seg_1.tau = (1.f - n_1C) * cot_th_1C;
        t.m_seg_2.tau = (1.f - n_2C) * cot_th_2C;

        // rho_theta
        if (c_perp != 0.f) {
            t.m_seg_1.rho_theta = -t.m_seg_1.tau / c_3D_1C;
            t.m_seg_2.rho_theta = -t.m_seg_2.tau / c_3D_2C;
        } else {
            t.m_seg_1.rho_theta = 0.f;
            t.m_seg_2.rho_theta = 0.f;
        }

        // theta
        t.m_seg_1.theta = theta_1C;
        t.m_seg_2.theta = theta_2C;

        // phi
        t.m_seg_1.phi = phi_01;
        t.m_seg_2.phi = phi_12;

        // nu
        t.m_seg_1.nu = (1.f - n_1C) * phi_1C;
        t.m_seg_2.nu = (1.f - n_2C) * phi_2C;

        // Triplet parameters

        t.m_phi_0 = 0.5f * (phi_1C * n_1C + phi_2C * n_2C);
        t.m_theta_0 = theta_2C - theta_1C +
                      ((1.f - n_2C) * cot_th_2C - (1.f - n_1C) * cot_th_1C);
        t.m_rho_phi = -0.5f * (t.m_seg_1.rho_phi + t.m_seg_2.rho_phi);
        t.m_rho_theta = t.m_seg_2.rho_theta - t.m_seg_1.rho_theta;
    }

    /// Helper function - Hit gradients
    ///
    /// Calculation of directional derivatives of
    /// triplet and segment parameters w.r.t hit
    /// position shifts
    ///
    /// @param t Triplet
    ///
    TRACCC_HOST_DEVICE void calculate_pos_derivs(
        triplet& t, edm::measurement_collection::const_device measurements) {

        // Hits shifted by multiplier * sigma in every direction
        scalar multiplier = 1.f;

        // Curvature of reference solution in bending plane
        scalar c_perp = -t.m_phi_0 / t.m_rho_phi;

        // Reserve space for derivative containers
        t.m_hd_phi.reserve(3u * m_max_dims);
        t.m_hd_theta.reserve(3u * m_max_dims);
        t.m_fseg_1.reserve(3u * m_max_dims);
        t.m_fseg_2.reserve(3u * m_max_dims);

        // Loop over hits
        for (unsigned hit = 0; hit < 3; ++hit) {

            // Over dimensions
            for (unsigned dir = 0; dir < m_max_dims; ++dir) {

                // Default derivative 0 for dimensions
                // which don't exist for this measurement
                if (dir >= measurements.at(t.m_meas_idx[hit]).dimensions()) {
                    t.m_hd_phi.push_back(0.f);
                    t.m_hd_theta.push_back(0.f);
                    continue;
                }

                // Make a copy to shift
                triplet t_shift{t};

                t_shift.shiftHit(hit, dir, multiplier, measurements,
                                 m_detector);

                // Calculate triplet and segment parameters
                linearize_triplet(t_shift, measurements);

                // Calculate derivatives

                // Triplet parameters
                scalar dphi = t_shift.m_phi_0 + c_perp * t_shift.m_rho_phi;
                scalar dtheta = (t_shift.m_theta_0 - t.m_theta_0) +
                                c_perp * (t_shift.m_rho_theta - t.m_rho_theta);
                t.m_hd_phi.push_back(dphi / multiplier);
                t.m_hd_theta.push_back(dtheta / multiplier);

                // Segment parameters
                t.m_fseg_1.emplace_back(
                    t_shift.m_seg_1.phi - t.m_seg_1.phi,
                    t_shift.m_seg_1.theta - t.m_seg_1.theta,
                    t_shift.m_seg_1.rho_phi - t.m_seg_1.rho_phi,
                    t_shift.m_seg_1.rho_theta - t.m_seg_1.rho_theta,
                    t_shift.m_seg_1.nu - t.m_seg_1.nu,
                    t_shift.m_seg_1.tau - t.m_seg_1.tau);

                t.m_fseg_2.emplace_back(
                    t_shift.m_seg_2.phi - t.m_seg_2.phi,
                    t_shift.m_seg_2.theta - t.m_seg_2.theta,
                    t_shift.m_seg_2.rho_phi - t.m_seg_2.rho_phi,
                    t_shift.m_seg_2.rho_theta - t.m_seg_2.rho_theta,
                    t_shift.m_seg_2.nu - t.m_seg_2.nu,
                    t_shift.m_seg_2.tau - t.m_seg_2.tau);
            }
        }
    }

    /// Helper function - Global Fit
    ///
    /// Global fit of processed hit triplets on track
    ///
    /// @param track track object being fitted
    /// @return fitted state at first measurement
    ///
    TRACCC_HOST_DEVICE
    typename edm::track_state_collection<algebra_type>::host::object_type
    do_global_fit(
        typename edm::track_collection<algebra_type>::host::proxy_type& track,
        edm::measurement_collection::const_device measurements) {

        // Allocate matrices with max possible sizes
        constexpr size_t max_nhits =
            20u;  // Assumption about max number of hits in track candidate
        constexpr size_t max_ntrips = max_nhits - 2u;
        constexpr size_t max_ndirs = m_max_dims * max_nhits;

        // Actual number in this track
        const size_t N_triplets = m_triplets.size();
        // assert(m_track_states.size() <= max_nhits);
        // assert(N_triplets == m_track_states.size() - 2u);
        assert(N_triplets + 2u <= max_nhits);

        // Make matrices/vectors

        // Triplet parameter vectors
        matrix_type<2u * max_ntrips, 1u> rho =
            matrix::zero<matrix_type<2u * max_ntrips, 1u>>();
        matrix_type<2u * max_ntrips, 1u> psi =
            matrix::zero<matrix_type<2u * max_ntrips, 1u>>();

        // Scattering covariance matrix
        matrix_type<2u * max_ntrips, 2u * max_ntrips> D_MS_inv =
            matrix::identity<matrix_type<2u * max_ntrips, 2u * max_ntrips>>();

        // Scaled hit gradient (Jacobian) matrix
        matrix_type<2u * max_ntrips, max_ndirs> Hd =
            matrix::zero<matrix_type<2u * max_ntrips, max_ndirs>>();

        // Fill matrices/vectors

        for (size_t i = 0; i < N_triplets; ++i) {

            const triplet& t_i = m_triplets[i];

            getter::element(rho, i, 0u) = t_i.m_rho_theta;
            getter::element(rho, i + max_ntrips, 0u) = t_i.m_rho_phi;

            getter::element(psi, i, 0u) = t_i.m_theta_0;
            getter::element(psi, i + max_ntrips, 0u) = t_i.m_phi_0;

            // Only update elements when linearization
            // has been done for this triplet
            // (reject on default value)
            if (t_i.m_sigma_MS != 0.f) {
                scalar sigma2_MS = t_i.m_sigma_MS * t_i.m_sigma_MS;
                scalar sin2_theta = math::sin(t_i.m_theta);
                sin2_theta *= sin2_theta;
                getter::element(D_MS_inv, i, i) = sigma2_MS;
                getter::element(D_MS_inv, i + max_ntrips, i + max_ntrips) =
                    sigma2_MS / sin2_theta;
            }

            // Hd_theta & phi
            // Loop unrolled over hits and uncertainty directions (assuming 2)
            // -----------------------------------------------
            // hd_theta and phi has the structure
            // hit0-dir0, hit0-dir1, hit1-dir0, hit1-dir1, ...
            // -----------------------------------------------

            // Theta derivatives
            // 1st Hit in triplet
            getter::element(Hd, i, i) = t_i.m_hd_theta[0u];
            getter::element(Hd, i, max_nhits + i) = t_i.m_hd_theta[1u];

            // 2nd Hit
            getter::element(Hd, i, i + 1u) = t_i.m_hd_theta[m_max_dims * 1u];
            getter::element(Hd, i, max_nhits + i + 1u) =
                t_i.m_hd_theta[m_max_dims * 1u + 1u];

            // 3rd Hit
            getter::element(Hd, i, i + 2u) = t_i.m_hd_theta[m_max_dims * 2u];
            getter::element(Hd, i, max_nhits + i + 2u) =
                t_i.m_hd_theta[m_max_dims * 2u + 1u];

            // Phi derivatives
            // 1st Hit
            getter::element(Hd, i + max_ntrips, i) = t_i.m_hd_phi[0u];
            getter::element(Hd, i + max_ntrips, max_nhits + i) =
                t_i.m_hd_phi[1u];

            // 2nd Hit
            getter::element(Hd, i + max_ntrips, i + 1u) =
                t_i.m_hd_phi[m_max_dims * 1u];
            getter::element(Hd, i + max_ntrips, max_nhits + i + 1u) =
                t_i.m_hd_phi[m_max_dims * 1u + 1u];

            // 3rd Hit
            getter::element(Hd, i + max_ntrips, i + 2u) =
                t_i.m_hd_phi[m_max_dims * 2u];
            getter::element(Hd, i + max_ntrips, max_nhits + i + 2u) =
                t_i.m_hd_phi[m_max_dims * 2u + 1u];

        }  // done filling

        // Triplet precision matrix
        // Note: diagonal elements in K_inv are 1
        // corresponding to unused 'objects', since
        // the same is true in D_MS_inv and those in
        // Hd * Hd^T are 0

        matrix_type<max_ndirs, 2u * max_ntrips> HdT = matrix::transpose(Hd);

        matrix_type<2u * max_ntrips, 2u * max_ntrips> K_inv =
            D_MS_inv + Hd * HdT;

        // Matrix inversion
        matrix_type<2u * max_ntrips, 2u * max_ntrips> K =
            matrix::inverse(K_inv);

        matrix_type<2u * max_ntrips, 1u> K_psi = K * psi;
        matrix_type<2u * max_ntrips, 1u> K_rho = K * rho;

        matrix_type<1u, 1u> rhoT_K_psi = matrix::transpose(rho) * K_psi;
        matrix_type<1u, 1u> rhoT_K_rho = matrix::transpose(rho) * K_rho;
        matrix_type<1u, 1u> psiT_K_psi = matrix::transpose(psi) * K_psi;

        // -------------------------------------------------- //
        // Calculation of curvature, uncertainty, fit quality //
        // -------------------------------------------------- //

        scalar c_3D = -1.f * getter::element(rhoT_K_psi, 0u, 0u) /
                      getter::element(rhoT_K_rho, 0u, 0u);

        scalar var_c_3D = 1.f / getter::element(rhoT_K_rho, 0u, 0u);

        scalar chi2 = getter::element(psiT_K_psi, 0u, 0u) +
                      c_3D * getter::element(rhoT_K_psi, 0u, 0u);

        // --------------------------------------------------- //
        // Calculation of correlation matrices & hit residuals //
        // --------------------------------------------------- //

        matrix_type<max_ndirs, 1u> HdT_K_psi = (HdT * K_psi);
        matrix_type<max_ndirs, 1u> HdT_K_rho = (HdT * K_rho);
        matrix_type<max_ndirs, 1u> HdT_K_rho_norm =
            (1.f / getter::element(rhoT_K_rho, 0u, 0u)) * HdT_K_rho;

        // Hit pull vector (local), all hits
        // Structure:
        // (h0,d0), (h1,d0), (h2,d0), ... , (hmax_nhits-1,d0),
        // (h0,d1), (h1,d1), (h2,d1), ... , (hmax_nhits-1,d1),
        // .
        // .
        // .
        // (h0,m_max_dims - 1), (h1,m_max_dims - 1), (h2,m_max_dims - 1), ... ,
        // (hmax_nhits-1,m_max_dims - 1)

        matrix_type<max_ndirs, 1u> hit_loc_pulls =
            HdT_K_psi -
            HdT_K_rho_norm *
                getter::element(rhoT_K_psi, 0u, 0u);  // to be stored ...

        // Hit matrix
        matrix_type<max_ndirs, max_ndirs> HdT_K_Hd =
            matrix::transpose(Hd) * K * Hd;

        // Hit correlation matrix (local)
        matrix_type<max_ndirs, max_ndirs> hit_loc_corrmatrix =
            matrix::identity<matrix_type<max_ndirs, max_ndirs>>() - HdT_K_Hd +
            HdT_K_rho_norm * matrix::transpose(HdT_K_rho);  // to be stored ...

        // Mixed hit position - curvature correlations (local)
        matrix_type<max_ndirs, 1u> hit_loc_corr_mixed =
            HdT_K_rho_norm;  // to be stored ...

        // ------------------------------------------------------- //
        // Calculation of track state vector and covariance matrix //
        // ------------------------------------------------------- //

        // At the first measurement surface (first triplet, first segement)

        const triplet& triplet_first = m_triplets[0u];

        scalar theta = triplet_first.m_seg_1.theta +
                       (triplet_first.m_seg_1.tau +
                        c_3D * triplet_first.m_seg_1.rho_theta);

        scalar phi = triplet_first.m_seg_1.phi -
                     0.5f * (triplet_first.m_seg_1.nu +
                             c_3D * triplet_first.m_seg_1.rho_phi);

        // grads_theta & phi have the same structure as
        // position and segment derivative containers
        constexpr unsigned grads_size = 3u * m_max_dims;
        std::array<scalar, grads_size> grads_theta{};
        std::array<scalar, grads_size> grads_phi{};

        // delta-theta & phi, delta-rho_theta & phi from hit position shifts
        scalar d_theta{};
        scalar d_phi{};
        scalar d_rho_theta{};
        scalar d_rho_phi{};

        for (unsigned i = 0; i < grads_size; ++i) {
            const segment_grad& grad = triplet_first.m_fseg_1[i];

            grads_theta[i] =
                grad.f_theta + grad.f_tau + c_3D * grad.f_rho_theta;
            grads_phi[i] =
                grad.f_phi + 0.5f * (grad.f_nu + c_3D * grad.f_rho_phi);

            unsigned idx = i / static_cast<unsigned>(m_max_dims) +
                           static_cast<unsigned>(max_nhits) *
                               (i % static_cast<unsigned>(m_max_dims));

            d_theta -= getter::element(hit_loc_pulls, idx, 0u) * grads_theta[i];
            d_phi -= getter::element(hit_loc_pulls, idx, 0u) * grads_phi[i];

            d_rho_theta +=
                getter::element(hit_loc_pulls, idx, 0u) * grad.f_rho_theta;
            d_rho_phi +=
                getter::element(hit_loc_pulls, idx, 0u) * grad.f_rho_phi;
        }

        // Apply track angle corrections from hit shifts
        theta += d_theta;
        phi = delta_phi(phi + d_phi, 0.f);

        // To rho_theta/phi
        scalar rho_theta = triplet_first.m_seg_1.rho_theta - d_rho_theta;
        scalar rho_phi = -0.5f * (triplet_first.m_seg_1.rho_phi - d_rho_phi);

        // Track parameters at first measurement surface
        detray::bound_parameters_vector<algebra_type> fitted_params;

        fitted_params.set_theta(theta);
        fitted_params.set_phi(phi);
        fitted_params.set_time(0.f);

        // Set momentum (q/p)

        const auto B_field = m_field.at(triplet_first.m_hit_pos[0][0u],
                                        triplet_first.m_hit_pos[0][1u],
                                        triplet_first.m_hit_pos[0][2u]);
        vector3 B_vec;
        B_vec[0u] = B_field[0u];
        B_vec[1u] = B_field[1u];
        B_vec[2u] = B_field[2u];

        // Units (in expression): B [T], p [MeV], c_3D [mm]
        // Expression: p = mom_conv * B / c_3D
        //-------------------------------------------------
        // Magentic field to be converted from traccc native to Tesla
        scalar p = mom_conv * detray::algebra::array::norm(B_vec) /
                   (c_3D * unit<scalar>::mm * unit<scalar>::T *
                    1000.f);  // GeV (native)

        const scalar q = 1.f;
        fitted_params.set_qop(q / p);

        // Post-fit position of first hit

        const auto& m0 = measurements.at(triplet_first.m_meas_idx[0]);

        point2 loc0{m0.local_position()[0], m0.local_position()[1]};
        point2 loc0_post_fit =
            loc0 + point2{getter::element(hit_loc_pulls, 0u, 0u),
                          getter::element(hit_loc_pulls, max_nhits,
                                          0u)};  // check sign: - ??

        fitted_params.set_bound_local(loc0_post_fit);

        // Track angle covariances
        scalar sum_thth{};
        scalar sum_phph{};
        scalar sum_thph{};
        // Track angle - curvature covariances
        scalar sum_c3D_th{};
        scalar sum_c3D_ph{};

        for (unsigned h1 = 0; h1 < grads_size; ++h1) {

            unsigned idx1 = h1 / static_cast<unsigned>(m_max_dims) +
                            static_cast<unsigned>(max_nhits) *
                                (h1 % static_cast<unsigned>(m_max_dims));

            sum_c3D_th +=
                getter::element(hit_loc_corr_mixed, idx1, 0u) * grads_theta[h1];
            sum_c3D_ph +=
                getter::element(hit_loc_corr_mixed, idx1, 0u) * grads_phi[h1];

            for (unsigned h2 = 0u; h2 < grads_size; ++h2) {

                unsigned idx2 = h2 / static_cast<unsigned>(m_max_dims) +
                                static_cast<unsigned>(max_nhits) *
                                    (h2 % static_cast<unsigned>(m_max_dims));

                sum_thth += getter::element(hit_loc_corrmatrix, idx1, idx2) *
                            grads_theta[h1] * grads_theta[h2];
                sum_thph += getter::element(hit_loc_corrmatrix, idx1, idx2) *
                            grads_theta[h1] * grads_phi[h2];
                sum_phph += getter::element(hit_loc_corrmatrix, idx1, idx2) *
                            grads_phi[h1] * grads_phi[h2];
            }
        }

        // Theta and Phi variance
        scalar var_theta = sum_thth - 2.f * sum_c3D_th * rho_theta +
                           var_c_3D * rho_theta * rho_theta;
        scalar var_phi = sum_phph - 2.f * sum_c3D_ph * rho_phi +
                         var_c_3D * rho_phi * rho_phi;

        // Theta-Phi covariance
        // scalar cov_theta_phi = sum_thph - sum_c3D_th * rho_phi - sum_c3D_ph *
        // rho_theta + rho_theta * rho_phi * var_c_3D;

        // Curvature-Theta/Phi covariance
        // scalar cov_c3D_theta = var_c_3D * rho_theta - sum_c3D_th;
        // scalar cov_c3D_phi = var_c_3D * rho_phi - sum_c3D_ph;

        // Covariance delta with theta/phi

        // variance of hit positions (for 1st hit)
        std::array<scalar, m_max_dims> var_pos_postFit;

        for (unsigned dir_i = 0; dir_i < m_max_dims; ++dir_i) {
            unsigned idx = dir_i * max_nhits;  // hit_i = 0
            var_pos_postFit[dir_i] =
                getter::element(hit_loc_corrmatrix, idx, idx) *
                measurements.at(triplet_first.m_meas_idx[0])
                    .local_variance()[dir_i];
        }

        // Covariance matrix
        matrix_type<6u, 6u> fit_cov{};

        getter::element(fit_cov, e_bound_loc0, e_bound_loc0) =
            var_pos_postFit[0u];
        getter::element(fit_cov, e_bound_loc1, e_bound_loc1) =
            var_pos_postFit[1u];

        getter::element(fit_cov, e_bound_phi, e_bound_phi) = var_phi;
        getter::element(fit_cov, e_bound_theta, e_bound_theta) = var_theta;

        getter::element(fit_cov, e_bound_qoverp, e_bound_qoverp) =
            var_c_3D * unit<scalar>::mm2 /
            math::pow((mom_conv * detray::algebra::array::norm(B_vec)) /
                          unit<scalar>::T,
                      2.f);
        getter::element(fit_cov, e_bound_time, e_bound_time) = 0;

        // Store final results
        track.chi2() = chi2;
        track.params().set_vector(fitted_params.vector());
        track.params().set_covariance(fit_cov);
        track.ndf() = m_meas_sum_dims - 5.f;
        track.pval() = prob(track.chi2(), track.ndf());

        if (chi2 > 0.f) {
            track.fit_outcome() = track_fit_outcome::SUCCESS;
        }

        // Only the smoothed parameters
        // at the first measurement are
        // used for performance plots
        // see fitting_performance_writer
        typename edm::track_state_collection<algebra_type>::host::object_type
            out_state{};

        out_state.smoothed_chi2() = chi2;
        out_state.smoothed_params().set_vector(fitted_params.vector());
        out_state.smoothed_params().set_covariance(fit_cov);
        out_state.set_hole(false);
        out_state.set_smoothed(true);

        return out_state;
    }

    /// Run the fitter - main fitting function
    ///
    /// @param track track object being fitted
    /// @return fitted track state at first surface
    ///
    TRACCC_HOST_DEVICE
    typename edm::track_state_collection<algebra_type>::host::object_type fit(
        typename edm::track_collection<algebra_type>::host::proxy_type& track,
        edm::measurement_collection::const_device measurements) {

        for (triplet& t : m_triplets) {

            linearize_triplet(t, measurements);
            calculate_pos_derivs(t, measurements);
        }

        auto first_state = do_global_fit(track, measurements);

        return first_state;
    }

    private:
    // Hard-coded material
    // TODO: get from surface after re-mapping
    scalar mat_scatter = 0.01f;

    // Conversion factor c_3D - p
    static constexpr scalar mom_conv = 0.299792458f;

    // Maximum number of uncertainty directions / hit
    static constexpr size_t m_max_dims = 2u;

    // Sum of measurement diensions
    scalar m_meas_sum_dims{};

    // Detector context type
    using context = typename detector_t::geometry_context;

    // Detector object
    const detector_t& m_detector;

    // Field object
    const bfield_t m_field;

    // Vector of triplets
    vecmem::vector<triplet> m_triplets;

    // Configuration object
    config_type m_cfg;
};

}  // namespace traccc

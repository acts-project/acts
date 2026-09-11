/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "traccc/definitions/primitives.hpp"
#include "traccc/definitions/qualifiers.hpp"
#include "traccc/definitions/track_parametrization.hpp"
#include "traccc/edm/measurement_helpers.hpp"
#include "traccc/finding/measurement_selector.hpp"
#include "traccc/fitting/details/regularize_covariance.hpp"
#include "traccc/fitting/status_codes.hpp"
#include "traccc/utils/logging.hpp"
#include "traccc/utils/matrix_helpers.hpp"
#include "traccc/utils/subspace.hpp"

namespace traccc {

/// Type unrolling functor for Kalman updating
template <typename algebra_t>
struct gain_matrix_updater {
  // Type declarations
  using size_type = detray::dindex_type<algebra_t>;
  template <size_type ROWS, size_type COLS>
  using matrix_type = detray::dmatrix<algebra_t, ROWS, COLS>;
  using bound_vector_type = traccc::bound_vector<algebra_t>;
  using bound_matrix_type = traccc::bound_matrix<algebra_t>;

  /// Gain matrix updater operation
  ///
  /// @brief Based on "Application of Kalman filtering to track and vertex
  /// fitting", R.Frühwirth, NIM A
  ///
  /// @param mask_group mask group that contains the mask of surface
  /// @param index mask index of surface
  /// @param trk_state track state of the surface
  /// @param bound_params bound parameter
  ///
  /// @return true if the update succeeds
  template <typename track_state_backend_t, typename measurement_backend_t>
  [[nodiscard]] TRACCC_HOST_DEVICE inline kalman_fitter_status operator()(
      typename edm::track_state<track_state_backend_t>& trk_state,
      const edm::measurement<measurement_backend_t>& meas,
      const bound_track_parameters<algebra_t>& bound_params,
      const measurement_selector::config& calib_cfg, const bool is_line) const {
    return this->operator()(trk_state.filtered_params(), meas, bound_params,
                            calib_cfg, is_line);
  }

  /// Gain matrix updater operation
  ///
  /// @brief Based on "Application of Kalman filtering to track and vertex
  /// fitting", R.Frühwirth, NIM A
  ///
  /// @param[out] filtered_params the filtered track vector and covariance
  /// @param[in] measurement the new measurement for the update
  /// @param[in] bound_params the predicted track parameters
  /// @param[in] is_line whether the measurement is on a line shaped surface
  ///
  /// @return kalman fitter status
  template <typename measurement_backend_t>
  [[nodiscard]] TRACCC_HOST_DEVICE inline kalman_fitter_status operator()(
      bound_track_parameters<algebra_t>& filtered_params,
      const edm::measurement<measurement_backend_t>& meas,
      const bound_track_parameters<algebra_t>& bound_params,
      const measurement_selector::config& calib_cfg, const bool is_line) const {
    static constexpr unsigned int D = 2;

    const unsigned int dim{meas.dimensions()};

    TRACCC_VERBOSE_HOST_DEVICE("Perform Kalman filtering...");
    TRACCC_VERBOSE_HOST_DEVICE("-> Measurement dim: %d", dim);

    assert(dim == 1u || dim == 2u);

    assert(!bound_params.is_invalid());
    assert(!bound_params.surface_link().is_invalid());

    TRACCC_DEBUG_HOST("Predicted param.: " << bound_params);

    // Predicted vector and covariance of bound track parameters
    const bound_vector_type& predicted_vec = bound_params.vector();
    const bound_matrix_type& predicted_cov = bound_params.covariance();

    // Measurement data on surface
    const matrix_type<D, 1> meas_local =
        measurement_selector::calibrated_measurement_position<algebra_t, D>(
            meas, calib_cfg);

    // Spatial resolution (Measurement covariance)
    const matrix_type<D, D> V =
        measurement_selector::calibrated_measurement_covariance<algebra_t, D>(
            meas, calib_cfg);

    const matrix_type<D, e_bound_size> H =
        measurement_selector::observation_model<algebra_t, D>(
            meas, bound_params, is_line);

    TRACCC_DEBUG_HOST("-> Predicted residual:\n"
                      << meas_local - H * predicted_vec);

    const matrix_type<e_bound_size, D> projected_cov =
        matrix::transposed_product<false, true>(predicted_cov, H);

    const matrix_type<D, D> M_inv =
        masked_inverse<algebra_t>(H * projected_cov + V, dim);

    // Kalman gain matrix
    const matrix_type<6, D> K = projected_cov * M_inv;

    TRACCC_DEBUG_HOST("-> H:\n" << H);
    TRACCC_DEBUG_HOST("-> K:\n" << K);

    // Calculate the filtered track parameters
    matrix_type<6, 1> filtered_vec =
        predicted_vec + K * (meas_local - H * predicted_vec);

    TRACCC_DEBUG_HOST("-> Filtered param:\n" << filtered_vec);

    // Return false if track is parallel to z-axis or phi is not finite
    if (!std::isfinite(getter::element(filtered_vec, e_bound_theta, 0))) {
      TRACCC_ERROR_HOST_DEVICE(
          "Theta is infinite after filtering (Matrix inversion)");
      return kalman_fitter_status::ERROR_INVERSION;
    }

    if (!std::isfinite(getter::element(filtered_vec, e_bound_phi, 0))) {
      TRACCC_ERROR_HOST_DEVICE(
          "Phi is infinite after filtering (Matrix inversion)");
      return kalman_fitter_status::ERROR_INVERSION;
    }

    if (math::fabs(getter::element(filtered_vec, e_bound_qoverp, 0)) == 0.f) {
      TRACCC_ERROR_HOST_DEVICE("q/p is zero after filtering");
      return kalman_fitter_status::ERROR_QOP_ZERO;
    }

    // Wrap the phi and theta angles in their valid ranges
    if (!normalize_angles<algebra_t>(filtered_vec)) {
      TRACCC_ERROR_HOST_DEVICE("Hit theta pole in filtering!");
      return kalman_fitter_status::ERROR_THETA_POLE;
    }

    // Some identity matrices
    // @TODO: Make constexpr work
    const auto I66 = matrix::identity<bound_matrix_type>();

    const matrix_type<6, 6> i_minus_kh = I66 - K * H;
    matrix_type<6, 6> filtered_cov =
        i_minus_kh * predicted_cov * matrix::transpose(i_minus_kh) +
        K * V * matrix::transpose(K);

    TRACCC_DEBUG_HOST("-> Filtered cov:\n" << filtered_cov);

    // Check the covariance for consistency
    // @TODO: Need to understand why negative variance happens
    if (constexpr traccc::scalar min_var{-0.01f};
        !details::regularize_covariance<algebra_t>(filtered_cov, min_var)) {
      TRACCC_ERROR_HOST_DEVICE("Negative variance after filtering");
      return kalman_fitter_status::ERROR_UPDATER_INVALID_COVARIANCE;
    }

    // Set the chi2 for this track and measurement
    filtered_params.set_vector(filtered_vec);
    filtered_params.set_covariance(filtered_cov);

    assert(!filtered_params.is_invalid());

    return kalman_fitter_status::SUCCESS;
  }
};

}  // namespace traccc

/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2024-2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "traccc/definitions/qualifiers.hpp"
#include "traccc/definitions/track_parametrization.hpp"
#include "traccc/edm/measurement_collection.hpp"
#include "traccc/edm/measurement_helpers.hpp"
#include "traccc/edm/track_state_collection.hpp"
#include "traccc/finding/measurement_selector.hpp"
#include "traccc/fitting/details/regularize_covariance.hpp"
#include "traccc/fitting/status_codes.hpp"
#include "traccc/utils/logging.hpp"
#include "traccc/utils/matrix_helpers.hpp"

namespace traccc {

/// Type unrolling functor for two-filters smoother
template <typename algebra_t>
struct two_filters_smoother {
  // Type declarations
  using size_type = detray::dindex_type<algebra_t>;
  template <size_type ROWS, size_type COLS>
  using matrix_type = detray::dmatrix<algebra_t, ROWS, COLS>;

  /// Two-filters smoother operation
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
      const edm::measurement<measurement_backend_t>& measurement,
      bound_track_parameters<algebra_t>& bound_params,
      const measurement_selector::config& calib_cfg, const bool is_line) const {
    static constexpr unsigned int D = 2;

    const unsigned int dim{measurement.dimensions()};

    assert(dim == 1u || dim == 2u);

    assert(!bound_params.is_invalid());
    assert(!bound_params.surface_link().is_invalid());
    assert(!trk_state.filtered_params().surface_link().is_invalid());
    assert(trk_state.filtered_params().surface_link() ==
           bound_params.surface_link());

    // Do not smoothe if the forward pass produced an error
    if (trk_state.filtered_params().is_invalid()) {
      TRACCC_ERROR_HOST_DEVICE("Filtered track state invalid");
      TRACCC_ERROR_HOST(trk_state.filtered_params());
      return kalman_fitter_status::ERROR_UPDATER_SKIPPED_STATE;
    }

    // Predicted vector of bound track parameters
    const matrix_type<e_bound_size, 1>& predicted_vec = bound_params.vector();

    // Predicted covaraince of bound track parameters
    const matrix_type<e_bound_size, e_bound_size>& predicted_cov =
        bound_params.covariance();

    const matrix_type<e_bound_size, e_bound_size> predicted_cov_inv =
        matrix::inverse(predicted_cov);
    const matrix_type<e_bound_size, e_bound_size> filtered_cov_inv =
        matrix::inverse(trk_state.filtered_params().covariance());

    // Eq (3.38) of "Pattern Recognition, Tracking and Vertex
    // Reconstruction in Particle Detectors"
    const matrix_type<e_bound_size, e_bound_size> smoothed_cov_inv =
        predicted_cov_inv + filtered_cov_inv;

    assert(matrix::determinant(smoothed_cov_inv) != 0.f);
    matrix_type<e_bound_size, e_bound_size> smoothed_cov =
        matrix::inverse(smoothed_cov_inv);

    // Check the covariance for consistency
    // @TODO: Need to understand why negative variance happens
    if (constexpr traccc::scalar min_var{-0.01f};
        !details::regularize_covariance<algebra_t>(smoothed_cov, min_var)) {
      TRACCC_ERROR_HOST_DEVICE("Negative variance after smoothing");
      return kalman_fitter_status::ERROR_SMOOTHER_INVALID_COVARIANCE;
    }

    // Eq (3.38) of "Pattern Recognition, Tracking and Vertex
    // Reconstruction in Particle Detectors"
    matrix_type<e_bound_size, 1u> smoothed_vec =
        smoothed_cov *
        (filtered_cov_inv * trk_state.filtered_params().vector() +
         predicted_cov_inv * predicted_vec);

    // Return false if track is parallel to z-axis or phi is not finite
    if (!std::isfinite(getter::element(smoothed_vec, e_bound_theta, 0))) {
      TRACCC_ERROR_HOST_DEVICE(
          "Theta is infinite after smoothing (Matrix inversion)");
      return kalman_fitter_status::ERROR_INVERSION;
    }

    if (!std::isfinite(getter::element(smoothed_vec, e_bound_phi, 0))) {
      TRACCC_ERROR_HOST_DEVICE(
          "Phi is infinite after smoothing (Matrix inversion)");
      return kalman_fitter_status::ERROR_INVERSION;
    }

    if (math::fabs(getter::element(smoothed_vec, e_bound_qoverp, 0)) == 0.f) {
      TRACCC_ERROR_HOST_DEVICE("q/p is zero after smoothing");
      return kalman_fitter_status::ERROR_QOP_ZERO;
    }

    // Wrap the phi and theta angles in their valid ranges
    if (!normalize_angles<algebra_t>(smoothed_vec)) {
      TRACCC_ERROR_HOST_DEVICE("Hit theta pole after smoothing!");
      return kalman_fitter_status::ERROR_THETA_POLE;
    }

    // Measurement data on surface
    const matrix_type<D, 1> meas_local =
        measurement_selector::calibrated_measurement_position<algebra_t, D>(
            measurement, calib_cfg);

    // Spatial resolution (Measurement covariance)
    const matrix_type<D, D> V =
        measurement_selector::calibrated_measurement_covariance<algebra_t, D>(
            measurement, calib_cfg);

    matrix_type<D, e_bound_size> H =
        measurement_selector::observation_model<algebra_t, D>(
            measurement, bound_params, is_line);

    const matrix_type<D, 1> residual_smt = meas_local - H * smoothed_vec;

    TRACCC_DEBUG_HOST("Predicted residual: " << meas_local - H * predicted_vec);

    // Eq (3.39) of "Pattern Recognition, Tracking and Vertex
    // Reconstruction in Particle Detectors"
    const matrix_type<D, D> R_smt =
        V - H * matrix::transposed_product<false, true>(smoothed_cov, H);

    const matrix_type<D, D> R_smt_inv = masked_inverse<algebra_t>(R_smt, dim);

    // Eq (3.40) of "Pattern Recognition, Tracking and Vertex
    // Reconstruction in Particle Detectors"
    const matrix_type<1, 1> chi2_smt =
        matrix::transposed_product<true, false>(residual_smt, R_smt_inv) *
        residual_smt;

    const scalar chi2_smt_value{getter::element(chi2_smt, 0, 0)};

    TRACCC_VERBOSE_HOST("Smoothed residual: " << residual_smt);
    TRACCC_DEBUG_HOST("R_smt:\n" << R_smt);
    TRACCC_DEBUG_HOST("R_smt_inv:\n" << R_smt_inv);
    TRACCC_VERBOSE_HOST_DEVICE("Smoothed chi2: %f", chi2_smt_value);

    if (chi2_smt_value < 0.f) {
      TRACCC_ERROR_HOST_DEVICE("Smoothed chi2 negative: %f", chi2_smt_value);

      // @TODO: Need to understand why negative chi2 happens
      if (chi2_smt_value < -10.f) {
        return kalman_fitter_status::ERROR_SMOOTHER_CHI2_NEGATIVE;
      }
    }

    if (!std::isfinite(chi2_smt_value)) {
      TRACCC_ERROR_HOST_DEVICE("Smoothed chi2 infinite");
      return kalman_fitter_status::ERROR_SMOOTHER_CHI2_NOT_FINITE;
    }

    /*************************************
     *  Set backward filtered parameter
     *************************************/

    const auto I66 =
        matrix::identity<matrix_type<e_bound_size, e_bound_size>>();
    const auto I_m = matrix::identity<matrix_type<D, D>>();

    const matrix_type<e_bound_size, D> projected_cov =
        matrix::transposed_product<false, true>(predicted_cov, H);

    const matrix_type<D, D> M_inv =
        masked_inverse<algebra_t>(H * projected_cov + V, dim);

    // Kalman gain matrix
    const matrix_type<6, D> K = projected_cov * M_inv;

    TRACCC_DEBUG_HOST("H:\n" << H);
    TRACCC_DEBUG_HOST("K:\n" << K);

    // Calculate the filtered track parameters
    matrix_type<6, 1> filtered_vec =
        predicted_vec + K * (meas_local - H * predicted_vec);

    // Return false if track is parallel to z-axis or phi is not finite
    if (!std::isfinite(getter::element(filtered_vec, e_bound_theta, 0))) {
      TRACCC_ERROR_HOST_DEVICE(
          "Theta is infinite after filering in smoother (Matrix "
          "inversion)");
      return kalman_fitter_status::ERROR_INVERSION;
    }

    if (!std::isfinite(getter::element(filtered_vec, e_bound_phi, 0))) {
      TRACCC_ERROR_HOST_DEVICE(
          "Phi is infinite after filering in smoother (Matrix "
          "inversion)");
      return kalman_fitter_status::ERROR_INVERSION;
    }

    if (math::fabs(getter::element(filtered_vec, e_bound_qoverp, 0)) == 0.f) {
      TRACCC_ERROR_HOST_DEVICE("q/p is zero after filering in smoother");
      return kalman_fitter_status::ERROR_QOP_ZERO;
    }

    // Wrap the phi and theta angles in their valid ranges
    if (!normalize_angles<algebra_t>(filtered_vec)) {
      TRACCC_ERROR_HOST_DEVICE("Hit theta pole after filtering in smoother!");
      return kalman_fitter_status::ERROR_THETA_POLE;
    }

    const matrix_type<6, 6> i_minus_kh = I66 - K * H;
    matrix_type<6, 6> filtered_cov =
        i_minus_kh * predicted_cov * matrix::transpose(i_minus_kh) +
        K * V * matrix::transpose(K);

    // Check the covariance for consistency
    // @TODO: Need to understand why negative variance happens
    if (constexpr traccc::scalar min_var{-0.01f};
        !details::regularize_covariance<algebra_t>(filtered_cov, min_var)) {
      TRACCC_ERROR_HOST_DEVICE("Negative variance after filtering");
      return kalman_fitter_status::ERROR_SMOOTHER_INVALID_COVARIANCE;
    }

    // Residual between measurement and (projected) filtered vector
    const matrix_type<D, 1> residual = meas_local - H * filtered_vec;

    // Calculate backward chi2
    const matrix_type<D, D> R = (I_m - H * K) * V;

    const matrix_type<D, D> R_inv = masked_inverse<algebra_t>(R, dim);

    const matrix_type<1, 1> chi2 =
        matrix::transposed_product<true, false>(residual, R_inv) * residual;

    const scalar chi2_val{getter::element(chi2, 0, 0)};

    TRACCC_VERBOSE_HOST("Filtered residual: " << residual);
    TRACCC_DEBUG_HOST("R:\n" << R);
    TRACCC_DEBUG_HOST("R_inv:\n" << R_inv);
    TRACCC_VERBOSE_HOST_DEVICE("Filtered chi2: %f", chi2_val);

    if (chi2_val < 0.f) {
      TRACCC_ERROR_HOST_DEVICE("Filtered chi2 negative: %f", chi2_val);
      return kalman_fitter_status::ERROR_SMOOTHER_CHI2_NEGATIVE;
    }

    if (!std::isfinite(chi2_val)) {
      TRACCC_ERROR_HOST_DEVICE("Filtered chi2 infinite");
      return kalman_fitter_status::ERROR_SMOOTHER_CHI2_NOT_FINITE;
    }

    // Update the smoothed track parameters
    trk_state.smoothed_params().set_vector(smoothed_vec);
    trk_state.smoothed_params().set_covariance(smoothed_cov);
    trk_state.smoothed_chi2() = getter::element(chi2_smt, 0, 0);
    trk_state.backward_chi2() = chi2_val;
    trk_state.set_smoothed();

    // Update the filtered track parameters
    bound_params.set_vector(filtered_vec);
    bound_params.set_covariance(filtered_cov);

    assert(!trk_state.smoothed_params().is_invalid());
    assert(!bound_params.is_invalid());

    return kalman_fitter_status::SUCCESS;
  }
};

}  // namespace traccc

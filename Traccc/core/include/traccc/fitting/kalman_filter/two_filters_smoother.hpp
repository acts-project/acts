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
#include "traccc/utils/concepts.hpp"
#include "traccc/utils/logging.hpp"
#include "traccc/utils/matrix_helpers.hpp"

// Detray include(s).
#include <detray/algebra/common/known_substructure_matrix.hpp>

namespace traccc {

/// Type unrolling functor for two-filters smoother
template <typename algebra_t>
struct two_filters_smoother {
  using scalar_type = detray::dscalar<algebra_t>;

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

    using scalar_t = detray::dscalar<algebra_t>;
    using cov_type =
        detray::ksm::matrix<track_covariance_substructure<e_bound_size>,
                            scalar_t>;
    using track_vec_type =
        detray::ksm::matrix<typename detray::ksm::make_dense_substructure<
                                e_bound_size, 1u>::canonical_type,
                            scalar_t>;
    using identity_type_bound =
        detray::ksm::matrix<typename detray::ksm::make_identity_substructure<
                                e_bound_size>::canonical_type,
                            scalar_t>;

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
    const concepts::column_matrix<e_bound_size, scalar_type> auto
        predicted_vec = track_vec_type::template from_dense<algebra_t>(
            bound_params.vector());
    const concepts::symmetric_matrix<e_bound_size, scalar_type> auto
        predicted_cov =
            symmetric_from_dense<cov_type>(bound_params.covariance());

    const concepts::symmetric_matrix<e_bound_size, scalar_type> auto
        predicted_cov_inv = symmetric_from_dense<cov_type>(
            matrix::inverse(predicted_cov.template to_dense<algebra_t>()));
    const concepts::symmetric_matrix<e_bound_size, scalar_type> auto
        filtered_cov_inv = symmetric_from_dense<cov_type>(
            matrix::inverse(trk_state.filtered_params().covariance()));

    // Eq (3.38) of "Pattern Recognition, Tracking and Vertex
    // Reconstruction in Particle Detectors"
    const concepts::symmetric_matrix<e_bound_size, scalar_type> auto
        smoothed_cov_inv = predicted_cov_inv + filtered_cov_inv;

    assert(matrix::determinant(
               smoothed_cov_inv.template to_dense<algebra_t>()) != 0.f);
    concepts::symmetric_matrix<e_bound_size, scalar_type> auto smoothed_cov =
        symmetric_from_dense<cov_type>(
            matrix::inverse(smoothed_cov_inv.template to_dense<algebra_t>()));

    // Check the covariance for consistency
    // @TODO: Need to understand why negative variance happens
    if (constexpr traccc::scalar min_var{-0.01f};
        !details::regularize_covariance(smoothed_cov, min_var)) {
      TRACCC_ERROR_HOST_DEVICE("Negative variance after smoothing");
      return kalman_fitter_status::ERROR_SMOOTHER_INVALID_COVARIANCE;
    }

    // Eq (3.38) of "Pattern Recognition, Tracking and Vertex
    // Reconstruction in Particle Detectors"
    concepts::column_matrix<e_bound_size, scalar_type> auto smoothed_vec =
        smoothed_cov *
        (filtered_cov_inv * track_vec_type::template from_dense<algebra_t>(
                                trk_state.filtered_params().vector()) +
         predicted_cov_inv * predicted_vec);

    // Return false if track is parallel to z-axis or phi is not finite
    if (!std::isfinite(getter::element<e_bound_theta, 0>(smoothed_vec))) {
      TRACCC_ERROR_HOST_DEVICE(
          "Theta is infinite after smoothing (Matrix inversion)");
      return kalman_fitter_status::ERROR_INVERSION;
    }

    if (!std::isfinite(getter::element<e_bound_phi, 0>(smoothed_vec))) {
      TRACCC_ERROR_HOST_DEVICE(
          "Phi is infinite after smoothing (Matrix inversion)");
      return kalman_fitter_status::ERROR_INVERSION;
    }

    if (math::fabs(getter::element<e_bound_qoverp, 0>(smoothed_vec)) == 0.f) {
      TRACCC_ERROR_HOST_DEVICE("q/p is zero after smoothing");
      return kalman_fitter_status::ERROR_QOP_ZERO;
    }

    // Wrap the phi and theta angles in their valid ranges
    if (!normalize_angles(smoothed_vec)) {
      TRACCC_ERROR_HOST_DEVICE("Hit theta pole after smoothing!");
      return kalman_fitter_status::ERROR_THETA_POLE;
    }

    // Measurement data on surface
    const concepts::column_matrix<D, scalar_type> auto meas_local =
        measurement_selector::calibrated_measurement_position<algebra_t, D>(
            measurement, calib_cfg);

    // Spatial resolution (Measurement covariance)
    const concepts::diagonal_matrix<D, scalar_type> auto V =
        measurement_selector::calibrated_measurement_covariance<algebra_t, D>(
            measurement, calib_cfg);

    const concepts::sized_matrix<D, e_bound_size, scalar_type> auto H =
        measurement_selector::observation_model<algebra_t, D>(
            measurement, bound_params, is_line);

    const concepts::column_matrix<D, scalar_type> auto residual_smt =
        meas_local - H * smoothed_vec;

    TRACCC_DEBUG_HOST(
        "Predicted residual: " << (meas_local - H * predicted_vec));

    // Eq (3.39) of "Pattern Recognition, Tracking and Vertex
    // Reconstruction in Particle Detectors"
    const concepts::symmetric_matrix<D, scalar_type> auto R_smt =
        V - H.congruence(smoothed_cov);

    const concepts::symmetric_matrix<D, scalar_type> auto R_smt_inv =
        masked_inverse(R_smt, dim);

    // Eq (3.40) of "Pattern Recognition, Tracking and Vertex
    // Reconstruction in Particle Detectors"
    const concepts::sized_matrix<1, 1, scalar_type> auto chi2_smt =
        residual_smt.transpose().congruence(R_smt_inv);

    const scalar chi2_smt_value{getter::element<0, 0>(chi2_smt)};

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

    const concepts::sized_matrix<e_bound_size, D, scalar_type> auto
        projected_cov = predicted_cov * H.transpose();

    const concepts::symmetric_matrix<D, scalar_type> auto M_inv =
        masked_inverse(H.congruence(predicted_cov) + V, dim);

    // Kalman gain matrix
    const concepts::sized_matrix<e_bound_size, D, scalar_type> auto K =
        projected_cov * M_inv;

    TRACCC_DEBUG_HOST("H:\n" << H);
    TRACCC_DEBUG_HOST("K:\n" << K);

    // Calculate the filtered track parameters
    concepts::column_matrix<e_bound_size, scalar_type> auto filtered_vec =
        predicted_vec + K * (meas_local - H * predicted_vec);

    // Return false if track is parallel to z-axis or phi is not finite
    if (!std::isfinite(getter::element<e_bound_theta, 0>(filtered_vec))) {
      TRACCC_ERROR_HOST_DEVICE(
          "Theta is infinite after filering in smoother (Matrix "
          "inversion)");
      return kalman_fitter_status::ERROR_INVERSION;
    }

    if (!std::isfinite(getter::element<e_bound_phi, 0>(filtered_vec))) {
      TRACCC_ERROR_HOST_DEVICE(
          "Phi is infinite after filering in smoother (Matrix "
          "inversion)");
      return kalman_fitter_status::ERROR_INVERSION;
    }

    if (math::fabs(getter::element<e_bound_qoverp, 0>(filtered_vec)) == 0.f) {
      TRACCC_ERROR_HOST_DEVICE("q/p is zero after filering in smoother");
      return kalman_fitter_status::ERROR_QOP_ZERO;
    }

    // Wrap the phi and theta angles in their valid ranges
    if (!normalize_angles(filtered_vec)) {
      TRACCC_ERROR_HOST_DEVICE("Hit theta pole after filtering in smoother!");
      return kalman_fitter_status::ERROR_THETA_POLE;
    }

    const concepts::sized_matrix<e_bound_size, e_bound_size, scalar_type> auto
        i_minus_kh = identity_type_bound::identity() - K * H;
    concepts::symmetric_matrix<e_bound_size, scalar_type> auto filtered_cov =
        i_minus_kh.congruence(predicted_cov) + K.congruence(V);

    // Check the covariance for consistency
    // @TODO: Need to understand why negative variance happens
    if (constexpr traccc::scalar min_var{-0.01f};
        !details::regularize_covariance(filtered_cov, min_var)) {
      TRACCC_ERROR_HOST_DEVICE("Negative variance after filtering");
      return kalman_fitter_status::ERROR_SMOOTHER_INVALID_COVARIANCE;
    }

    // Residual between measurement and (projected) filtered vector
    const concepts::column_matrix<D, scalar_type> auto residual =
        meas_local - H * filtered_vec;

    // Calculate backward chi2
    //
    // R is the covariance of the filtered residual, and the usual way to
    // write it is (I - H*K)*V. That form hides a cancellation. Since
    //
    // K = C*H^T*M^-1 ; and
    // M = H*C*H^T + V
    //
    // It follows that:
    //
    // H*K = (H*C*H^T)*M^-1 = (M - V)*M^-1 = I - V*M^-1
    //
    // Such that:
    //
    // I - H*K = I - (I - V*M^-1) = V*M^-1
    //
    // And thus:
    //
    // (I - H*K)*V = (V*M^-1)*V
    //
    // We use the collapsed form for two reasons. It is a congruence, so it
    // is symmetric because of its shape. It also doesn't subtract anything,
    // nothing, so it is less susceptible to catastrophic cancellation.
    const concepts::symmetric_matrix<D, scalar_type> auto R =
        V.congruence(M_inv);

    const concepts::symmetric_matrix<D, scalar_type> auto R_inv =
        masked_inverse(R, dim);

    const concepts::sized_matrix<1, 1, scalar_type> auto chi2 =
        residual.transpose().congruence(R_inv);

    const scalar chi2_val{getter::element<0, 0>(chi2)};

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
    trk_state.smoothed_params().set_vector(
        smoothed_vec.template to_dense<algebra_t>());
    trk_state.smoothed_params().set_covariance(
        smoothed_cov.template to_dense<algebra_t>());
    trk_state.smoothed_chi2() = getter::element<0, 0>(chi2_smt);
    trk_state.backward_chi2() = chi2_val;
    trk_state.set_smoothed();

    // Update the filtered track parameters
    bound_params.set_vector(filtered_vec.template to_dense<algebra_t>());
    bound_params.set_covariance(filtered_cov.template to_dense<algebra_t>());

    assert(!trk_state.smoothed_params().is_invalid());
    assert(!bound_params.is_invalid());

    return kalman_fitter_status::SUCCESS;
  }
};

}  // namespace traccc

/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "traccc/definitions/primitives.hpp"
#include "traccc/definitions/qualifiers.hpp"
#include "traccc/definitions/track_parametrization.hpp"
#include "traccc/edm/measurement_collection.hpp"
#include "traccc/edm/measurement_helpers.hpp"
#include "traccc/edm/track_parameters.hpp"
#include "traccc/utils/logging.hpp"
#include "traccc/utils/matrix_helpers.hpp"
#include "traccc/utils/subspace.hpp"

// Detray include(s).
#include <detray/algebra/common/known_substructure_matrix.hpp>

// System include(s)
#include <limits>

namespace traccc {

/// Potential next measurement to be added to track
struct candidate_measurement {
  unsigned int meas_idx{std::numeric_limits<unsigned int>::max()};
  float chi2{std::numeric_limits<float>::max()};

  /// Define comparisons
  constexpr bool operator<=>(const candidate_measurement& other) const =
      default;
};

/// Associate a measurement to a candidate track
struct measurement_selector {
  template <detray::concepts::algebra A, unsigned int ROWS, unsigned int COLS>
  using matrix_t = detray::dmatrix<A, ROWS, COLS>;

  // Where to get the calibration from
  struct config { /*TODO: implement calibration handling*/
  };

  /// Get the observation model for a given measurement
  ///
  /// @param measurement the measurement
  /// @param bound_params predicted bound track parameters
  /// @param is_line whether the measurement belong to a line surface
  ///
  /// @returns the projection matrix H
  ///
  /// @note A measurement subspace only ever selects @c e_bound_loc0 and
  /// @c e_bound_loc1, so the result declares every further column a
  /// structural zero. Call @c to_dense on it where a dense matrix is needed.
  template <detray::concepts::algebra algebra_t, unsigned int D,
            typename measurement_backend_t>
  TRACCC_HOST_DEVICE static detray::ksm::matrix<
      observation_model_substructure<D, e_bound_size>,
      detray::dscalar<algebra_t>>
  observation_model(const edm::measurement<measurement_backend_t>& measurement,
                    const bound_track_parameters<algebra_t>& bound_params,
                    const bool is_line) {
    // Oservation model: Subspace of measurement space for this measurement
    subspace<algebra_t, e_bound_size> subs(measurement.subspace());

    // Flip the sign of projector matrix element in case the first element
    // of a line measurement is negative
    if (is_line && bound_params.bound_local()[e_bound_loc0] < 0) {
      subs.set_sign(0, true);
    }

    if (measurement.dimensions() == 1) {
      subs.set_invalid(1);
    }

    const auto H = subs.template projector<D, 2u>();

    TRACCC_DEBUG_HOST("--> Observation model (H):\n" << H);

    return H;
  }

  /// Get the calibrated measurement position
  ///
  /// @param measurement the measurement
  /// @param cfg how to apply calibrations
  ///
  /// @returns the calibrated measurement position
  ///
  /// @note The result is a known-substructure matrix. Call @c to_dense on it
  /// where a dense matrix is needed.
  template <detray::concepts::algebra algebra_t, unsigned int D,
            typename measurement_backend_t>
  TRACCC_HOST_DEVICE static detray::ksm::matrix<
      typename detray::ksm::make_dense_substructure<D, 1u>::canonical_type,
      detray::dscalar<algebra_t>>
  calibrated_measurement_position(
      const edm::measurement<measurement_backend_t>& measurement,
      const config& /*cfg*/) {
    // Measurement local position on surface. The EDM stores it densely, so
    // the one conversion lives here rather than at each of the callers.
    detray::dmatrix<algebra_t, D, 1> meas_local;
    edm::get_measurement_local<algebra_t>(measurement, meas_local);

    TRACCC_DEBUG_HOST("--> Measurement position (uncalibrated):\n"
                      << meas_local);

    assert((measurement.dimensions() > 1) ||
           (getter::element(meas_local, 1u, 0u) == 0.f));

    return detray::ksm::matrix<
        typename detray::ksm::make_dense_substructure<D, 1u>::canonical_type,
        detray::dscalar<algebra_t>>::template from_dense<algebra_t>(meas_local);
  }

  /// Get the calibrated measurement covariance
  ///
  /// @param measurement the measurement
  /// @param cfg how to apply calibrations
  ///
  /// @returns the calibrated measurement covariance
  ///
  /// @note The result declares itself diagonal, because
  /// @c edm::get_measurement_covariance writes the off-diagonal elements as
  /// exact zeros. Call @c to_dense on it where a dense matrix is needed.
  template <detray::concepts::algebra algebra_t, unsigned int D,
            typename measurement_backend_t>
  TRACCC_HOST_DEVICE static detray::ksm::matrix<
      measurement_covariance_substructure<D>, detray::dscalar<algebra_t>>
  calibrated_measurement_covariance(
      const edm::measurement<measurement_backend_t>& measurement,
      const config& /*cfg*/) {
    // Measurement covariance. The EDM stores it densely, so the one
    // conversion lives here rather than at each of the callers.
    detray::dmatrix<algebra_t, D, D> V;
    edm::get_measurement_covariance<algebra_t>(measurement, V);

    if (measurement.dimensions() == 1) {
      getter::element(V, 1u, 1u) =
          std::numeric_limits<detray::dscalar<algebra_t>>::max();
    }

    TRACCC_DEBUG_HOST("--> Measurement covariance (uncalibrated):\n" << V);

    return detray::ksm::matrix<
        measurement_covariance_substructure<D>,
        detray::dscalar<algebra_t>>::template from_dense<algebra_t>(V);
  }

  /// Calculate the predicted chi2
  ///
  /// @brief Based on "Application of Kalman filtering to track and vertex
  /// fitting", R.Frühwirth, NIM A
  ///
  /// @param measurement the measurement
  /// @param bound_params predicted bound track parameters
  /// @param cfg the calibration configuration
  /// @param is_line whether the measurement belong to a line surface
  ///
  /// @returns the predicted chi2 of the calibrated measurement
  template <typename measurement_backend_t, detray::concepts::algebra algebra_t>
  TRACCC_HOST_DEVICE static detray::dscalar<algebra_t> predicted_chi2(
      const edm::measurement<measurement_backend_t>& measurement,
      const bound_track_parameters<algebra_t>& bound_params, const config& cfg,
      const bool is_line) {
    using scalar_t = detray::dscalar<algebra_t>;

    // Measurement maximal dimension
    constexpr unsigned int D = 2;

    namespace ksm = detray::ksm;

    // The shapes this computation reads its inputs into. The observation
    // model needs none, because it already arrives structured.
    using cov_type =
        ksm::matrix<track_covariance_substructure<e_bound_size>, scalar_t>;
    using track_vec_type = ksm::matrix<
        typename ksm::make_dense_substructure<e_bound_size, 1u>::canonical_type,
        scalar_t>;

    TRACCC_VERBOSE_HOST_DEVICE("--> dim: %d", measurement.dimensions());

    assert(measurement.dimensions() == 1u || measurement.dimensions() == 2u);

    assert(!bound_params.is_invalid());
    assert(!bound_params.surface_link().is_invalid());

    // Get calibrated measurement and covariance
    const auto meas_local =
        calibrated_measurement_position<algebra_t, D>(measurement, cfg);

    const auto V =
        calibrated_measurement_covariance<algebra_t, D>(measurement, cfg);

    // Project the predicted covariance to the observation. Only the upper
    // triangle of the covariance is read, because the two mirrored elements
    // of a transported covariance are equal only up to rounding.
    const auto H =
        observation_model<algebra_t, D>(measurement, bound_params, is_line);

    const auto predicted_cov =
        symmetric_from_dense<cov_type>(bound_params.covariance());

    // H.congruence(C) is H*C*H^T, and it carries the symmetry of the result
    // in its type.
    const auto R = H.congruence(predicted_cov) + V;

    const auto R_inv = masked_inverse(R, measurement.dimensions());

    TRACCC_DEBUG_HOST("--> R:\n" << R);
    TRACCC_DEBUG_HOST("--> R_inv:\n" << R_inv);

    // Residual between measurement and (projected) vector (innovation)
    const auto residual =
        meas_local - H * track_vec_type::template from_dense<algebra_t>(
                             bound_params.vector());

    TRACCC_DEBUG_HOST("--> Predicted residual:\n" << residual);

    const auto pred_chi2 = (residual.transpose() * R_inv) * residual;

    const scalar_t pred_chi2_val{pred_chi2.template at<0, 0>()};

    if (!std::isfinite(pred_chi2_val)) {
      TRACCC_WARNING_HOST_DEVICE("Infinite predicted chi2 value!");
    } else if (pred_chi2_val < 0.f) {
      TRACCC_WARNING_HOST_DEVICE("Negative predicted chi2 value!");
    } else {
      TRACCC_VERBOSE_HOST_DEVICE("--> chi2: %.10e", pred_chi2_val);
    }

    return pred_chi2_val;
  }

  /// Measurement selection (optimal)
  ///
  /// @param bound_params predicted bound track parameters
  /// @param measurements the measurement container
  /// @param meas_range contains the index ranges into the measurements
  /// @param cfg the calibration configuration
  /// @param is_line whether the measurement belong to a line surface
  ///
  /// @returns the optimal candidate measurement for the input params
  template <detray::concepts::algebra algebra_t>
  TRACCC_HOST_DEVICE static candidate_measurement find_optimal_measurement(
      const bound_track_parameters<algebra_t>& bound_params,
      const typename edm::measurement_collection::const_device& measurements,
      vecmem::device_vector<unsigned int> meas_ranges, const config& cfg,
      const bool is_line) {
    using scalar_t = detray::dscalar<algebra_t>;

    // The optimal candidate
    candidate_measurement cand{};

    // Iterate over the measurements for this surface
    const unsigned int sf_idx{bound_params.surface_link().index()};
    const unsigned int lo{sf_idx == 0u ? 0u : meas_ranges[sf_idx - 1]};
    const unsigned int up{meas_ranges[sf_idx]};

    TRACCC_VERBOSE_HOST_DEVICE("Have %d measurement(s) on surface %d...",
                               up - lo, sf_idx);

    // Find the best fitting measurement by prediced chi2
    // TODO: Load balancing
    for (unsigned int meas_idx = lo; meas_idx < up; meas_idx++) {
      TRACCC_VERBOSE_HOST_DEVICE("-> measurement %d:", meas_idx);

      // Predicted chi2
      const scalar_t chi2 = measurement_selector::predicted_chi2(
          measurements.at(meas_idx), bound_params, cfg, is_line);

      // Check predicted chi2 cut
      if (chi2 < cand.chi2 && chi2 >= 0.f) {
        cand = {meas_idx, static_cast<float>(chi2)};
        // Found optimal
        if (cand.chi2 <= std::numeric_limits<scalar_t>::epsilon()) {
          return cand;
        }
      }
    }

    return cand;
  }

  /// Measurement selection (collection of compatible measurements)
  ///
  /// @param bound_params predicted bound track parameters
  /// @param measurements the measurement container
  /// @param meas_range contains the index ranges into the measurements
  /// @param cfg the calibration configuration
  /// @param is_line whether the measurement belong to a line surface
  ///
  /// @returns a collection of compatible measurements, sorted by pred.
  /// chi2
  template <detray::concepts::algebra algebra_t>
  TRACCC_HOST_DEVICE static vecmem::vector<candidate_measurement>
  find_compatible_measurements(
      const bound_track_parameters<algebra_t>& /*bound_params*/,
      const typename edm::measurement_collection::const_device&
      /*measurements*/,
      vecmem::device_vector<unsigned int> /*meas_ranges*/,
      const config& /*cfg*/, const bool /*is_line*/) {
    /* TODO: Implement*/
    assert(false);
    return {};
  }
};

}  // namespace traccc

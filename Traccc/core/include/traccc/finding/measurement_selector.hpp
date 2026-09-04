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
  template <detray::concepts::algebra algebra_t, unsigned int D,
            typename measurement_backend_t>
  TRACCC_HOST_DEVICE static detray::dmatrix<algebra_t, D, e_bound_size>
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

    detray::dmatrix<algebra_t, D, 1> meas_local;
    edm::get_measurement_local<algebra_t>(measurement, meas_local);
    if (measurement.dimensions() == 1) {
      subs.set_invalid(1);
    }

    const detray::dmatrix<algebra_t, D, e_bound_size> H =
        subs.template projector<D>();

    TRACCC_DEBUG_HOST("--> Observation model (H):\n" << H);

    return H;
  }

  /// Get the calibrated measurement position
  ///
  /// @param measurement the measurement
  /// @param cfg how to apply calibrations
  ///
  /// @returns the projection matrix H
  template <detray::concepts::algebra algebra_t, unsigned int D,
            typename measurement_backend_t>
  TRACCC_HOST_DEVICE static detray::dmatrix<algebra_t, D, 1>
  calibrated_measurement_position(
      const edm::measurement<measurement_backend_t>& measurement,
      const config& /*cfg*/) {
    // Measurement local position on surface
    detray::dmatrix<algebra_t, D, 1> meas_local;
    edm::get_measurement_local<algebra_t>(measurement, meas_local);

    TRACCC_DEBUG_HOST("--> Measurement position (uncalibrated):\n"
                      << meas_local);

    assert((measurement.dimensions() > 1) ||
           (getter::element(meas_local, 1u, 0u) == 0.f));

    return meas_local;
  }

  /// Get the calibrated measurement covariance
  ///
  /// @param measurement the measurement
  /// @param cfg how to apply calibrations
  ///
  /// @returns the projection matrix H
  template <detray::concepts::algebra algebra_t, unsigned int D,
            typename measurement_backend_t>
  TRACCC_HOST_DEVICE static detray::dmatrix<algebra_t, D, D>
  calibrated_measurement_covariance(
      const edm::measurement<measurement_backend_t>& measurement,
      const config& /*cfg*/) {
    // Measurement covariance
    detray::dmatrix<algebra_t, D, D> V;
    edm::get_measurement_covariance<algebra_t>(measurement, V);

    detray::dmatrix<algebra_t, D, 1> meas_local;
    edm::get_measurement_local<algebra_t>(measurement, meas_local);
    if (measurement.dimensions() == 1) {
      getter::element(V, 1u, 1u) =
          std::numeric_limits<detray::dscalar<algebra_t>>::max();
    }

    TRACCC_DEBUG_HOST("--> Measurement covariance (uncalibrated):\n" << V);

    return V;
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

    TRACCC_VERBOSE_HOST_DEVICE("--> dim: %d", measurement.dimensions());

    assert(measurement.dimensions() == 1u || measurement.dimensions() == 2u);

    assert(!bound_params.is_invalid());
    assert(!bound_params.surface_link().is_invalid());

    // Get calibrated measurement and covariance
    const matrix_t<algebra_t, D, 1> meas_local =
        calibrated_measurement_position<algebra_t, D>(measurement, cfg);

    const matrix_t<algebra_t, D, D> V =
        calibrated_measurement_covariance<algebra_t, D>(measurement, cfg);

    // Project the predicted covariance to the observation
    const matrix_t<algebra_t, D, e_bound_size> H =
        observation_model<algebra_t, D>(measurement, bound_params, is_line);

    const matrix_t<algebra_t, D, D> R =
        H * matrix::transposed_product<false, true>(bound_params.covariance(),
                                                    H) +
        V;

    const matrix_t<algebra_t, D, D> R_inv =
        masked_inverse<algebra_t>(R, measurement.dimensions());

    TRACCC_DEBUG_HOST("--> R:\n" << R);
    TRACCC_DEBUG_HOST("--> R_inv:\n" << R_inv);

    // Residual between measurement and (projected) vector (innovation)
    const matrix_t<algebra_t, D, 1> residual =
        meas_local - H * bound_params.vector();

    TRACCC_DEBUG_HOST("--> Predicted residual:\n" << residual);

    const matrix_t<algebra_t, 1, 1> pred_chi2 =
        matrix::transposed_product<true, false>(residual, R_inv) * residual;

    const scalar_t pred_chi2_val{getter::element(pred_chi2, 0, 0)};

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

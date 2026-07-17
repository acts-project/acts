/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2025-2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Local include(s).
#include "traccc/utils/detray_conversion.hpp"

namespace traccc::edm {

template <detray::concepts::algebra algebra_t, typename measurement_backend_t>
TRACCC_HOST_DEVICE detray::dpoint2D<algebra_t> get_measurement_local(
    const edm::measurement<measurement_backend_t>& meas) {

    return utils::to_dpoint2D<algebra_t>(meas.local_position());
}

template <detray::concepts::algebra algebra_t, typename measurement_backend_t,
          std::integral size_t, size_t D>
TRACCC_HOST_DEVICE void get_measurement_local(
    const edm::measurement<measurement_backend_t>& meas,
    detray::dmatrix<algebra_t, D, 1>& pos) {

    static_assert(((D == 1u) || (D == 2u)),
                  "The measurement dimension must be 1 or 2");

    assert((meas.subspace()[0] == detray::e_bound_loc0) ||
           (meas.subspace()[0] == detray::e_bound_loc1));

    switch (meas.subspace()[0]) {
        case detray::e_bound_loc0:
            getter::element(pos, 0, 0) =
                static_cast<typename algebra_t::value_type>(
                    meas.local_position()[0]);
            if constexpr (D == 2u) {
                getter::element(pos, 1, 0) =
                    static_cast<typename algebra_t::value_type>(
                        meas.local_position()[1]);
            }
            break;
        case detray::e_bound_loc1:
            getter::element(pos, 0, 0) =
                static_cast<typename algebra_t::value_type>(
                    meas.local_position()[1]);
            if constexpr (D == 2u) {
                getter::element(pos, 1, 0) =
                    static_cast<typename algebra_t::value_type>(
                        meas.local_position()[0]);
            }
            break;
        default:
#if defined(__GNUC__)
            __builtin_unreachable();
#endif
    }
}

template <detray::concepts::algebra algebra_t, typename measurement_backend_t>
TRACCC_HOST_DEVICE detray::dvector2D<algebra_t> get_measurement_variance(
    const edm::measurement<measurement_backend_t>& meas) {

    return utils::to_dpoint2D<algebra_t>(meas.local_variance());
}

template <detray::concepts::algebra algebra_t, typename measurement_backend_t,
          std::integral size_t, size_t D>
TRACCC_HOST_DEVICE void get_measurement_covariance(
    const edm::measurement<measurement_backend_t>& meas,
    detray::dmatrix<algebra_t, D, D>& cov) {

    static_assert(((D == 1u) || (D == 2u)),
                  "The measurement dimension must be 1 or 2");

    assert((meas.subspace()[0] == detray::e_bound_loc0) ||
           (meas.subspace()[0] == detray::e_bound_loc1));

    switch (meas.subspace()[0]) {
        case detray::e_bound_loc0:
            getter::element(cov, 0, 0) =
                static_cast<typename algebra_t::value_type>(
                    meas.local_variance()[0]);
            if constexpr (D == 2u) {
                getter::element(cov, 0, 1) = 0.f;
                getter::element(cov, 1, 0) = 0.f;
                getter::element(cov, 1, 1) =
                    static_cast<typename algebra_t::value_type>(
                        meas.local_variance()[1]);
            }
            break;
        case detray::e_bound_loc1:
            getter::element(cov, 0, 0) =
                static_cast<typename algebra_t::value_type>(
                    meas.local_variance()[1]);
            if constexpr (D == 2u) {
                getter::element(cov, 0, 1) = 0.f;
                getter::element(cov, 1, 0) = 0.f;
                getter::element(cov, 1, 1) =
                    static_cast<typename algebra_t::value_type>(
                        meas.local_variance()[0]);
            }
            break;
        default:
#if defined(__GNUC__)
            __builtin_unreachable();
#endif
    }
}

}  // namespace traccc::edm

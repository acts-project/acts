/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "traccc/definitions/qualifiers.hpp"

// Detray include(s).
#include <detray/definitions/algebra.hpp>

// System include(s).
#include <array>

namespace traccc::utils {

/// Convert an algebra specific 2D point to a float array of size 2
///
/// @tparam ALGEBRA_TYPE The algebra type of the input point
/// @param point The input 2D point to be converted
/// @return An @c std::array of size 2 containing the converted float values
///
template <detray::concepts::algebra ALGEBRA_TYPE>
TRACCC_HOST_DEVICE std::array<float, 2u> to_float_array(
    const detray::dpoint2D<ALGEBRA_TYPE>& point);

/// Convert an algebra specific 3D point to a float array of size 3
///
/// @tparam ALGEBRA_TYPE The algebra type of the input point
/// @param point The input 3D point to be converted
/// @return An @c std::array of size 3 containing the converted float values
///
template <detray::concepts::algebra ALGEBRA_TYPE>
TRACCC_HOST_DEVICE std::array<float, 3u> to_float_array(
    const detray::dpoint3D<ALGEBRA_TYPE>& point);

/// Convert a float array of size 2 to an algebra specific 2D point
///
/// @tparam ALGEBRA_TYPE The algebra type of the output point
/// @param arr The input @c std::array of size 2 containing the float values
/// @return A 2D point of type @c detray::dpoint2D with the converted values
///
template <detray::concepts::algebra ALGEBRA_TYPE>
TRACCC_HOST_DEVICE detray::dpoint2D<ALGEBRA_TYPE> to_dpoint2D(
    const std::array<float, 2u>& arr);

/// Convert a float array of size 3 to an algebra specific 3D point
///
/// @tparam ALGEBRA_TYPE The algebra type of the output point
/// @param arr The input @c std::array of size 3 containing the float values
/// @return A 3D point of type @c detray::dpoint3D with the converted values
///
template <detray::concepts::algebra ALGEBRA_TYPE>
TRACCC_HOST_DEVICE detray::dpoint3D<ALGEBRA_TYPE> to_dpoint3D(
    const std::array<float, 3u>& arr);

}  // namespace traccc::utils

// Implementation include(s).
#include "traccc/utils/impl/detray_conversion.ipp"

/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2025 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "traccc/definitions/common.hpp"
#include "traccc/definitions/math.hpp"
#include "traccc/definitions/primitives.hpp"
#include "traccc/definitions/qualifiers.hpp"

/// @see
/// https://github.com/acts-project/acts/blob/8098e6953ac35771c34a4e3b13dbfee50869a0c1/Core/include/Acts/Utilities/detail/periodic.hpp
namespace traccc::detail {

/// Wrap a periodic value back into the nominal range.
TRACCC_HOST_DEVICE
constexpr traccc::scalar wrap_periodic(traccc::scalar value,
                                       traccc::scalar start,
                                       traccc::scalar range) {
    // only wrap if really necessary
    const traccc::scalar diff{value - start};
    return ((0 <= diff) && (diff < range))
               ? value
               : (value - range * math::floor(diff / range));
}

/// Calculate the equivalent angle in the [-pi, pi) range.
TRACCC_HOST_DEVICE
constexpr traccc::scalar wrap_phi(traccc::scalar phi) {
    constexpr traccc::scalar PI{traccc::constant<traccc::scalar>::pi};
    constexpr traccc::scalar TWOPI{2.f * traccc::constant<traccc::scalar>::pi};
    return wrap_periodic(phi, -PI, TWOPI);
}

/// Calculate the equivalent angle in the [0, 2*pi) range.
TRACCC_HOST_DEVICE
constexpr traccc::scalar wrap_theta(traccc::scalar theta) {
    constexpr traccc::scalar TWOPI{2.f * traccc::constant<traccc::scalar>::pi};
    return wrap_periodic(theta, 0.f, TWOPI);
}

/// Ensure both phi and theta direction angles are within the allowed range.
///
/// @param[in] phi Transverse direction angle
/// @param[in] theta Longitudinal direction angle
/// @return pair<phi,theta> containing the updated angles
///
/// The phi angle is truly cyclic, i.e. all values outside the nominal range
/// [-pi,pi) have a corresponding value inside nominal range, independent from
/// the theta angle. The theta angle is more complicated. Imagine that the two
/// angles describe a position on the unit sphere. If theta moves outside its
/// nominal range [0,pi], we are moving over one of the two poles of the unit
/// sphere along the great circle defined by phi. The angles still describe a
/// valid position on the unit sphere, but to describe it with angles within
/// their nominal range, both phi and theta need to be updated; when moving over
/// the poles, phi needs to be flipped by 180degree to allow theta to remain
/// within its nominal range.
constexpr std::pair<traccc::scalar, traccc::scalar> wrap_phi_theta(
    traccc::scalar phi, traccc::scalar theta) {
    constexpr traccc::scalar PI{traccc::constant<traccc::scalar>::pi};
    constexpr traccc::scalar TWOPI{2.f * traccc::constant<traccc::scalar>::pi};

    // wrap to [0,2pi). while the nominal range of theta is [0,pi], it is
    // periodic, i.e. describes identical positions, in the full [0,2pi) range.
    // moving it first to the periodic range simplifies further steps as the
    // possible range of theta becomes fixed.
    theta = wrap_theta(theta);
    if (PI < theta) {
        // theta is in the second half of the great circle and outside its
        // nominal range. need to change both phi and theta to be within range.
        phi += PI;
        theta = TWOPI - theta;
    }

    return {wrap_phi(phi), theta};
}

}  // namespace traccc::detail

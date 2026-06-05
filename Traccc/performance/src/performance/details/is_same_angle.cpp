/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Library include(s).
#include "traccc/performance/details/is_same_angle.hpp"

// System include(s).
#include <cmath>
#include <numbers>

namespace {

inline traccc::scalar wrap_to_pi(traccc::scalar phi) {

    // Make sure that we only use the precision necessary.
    static constexpr traccc::scalar PI = std::numbers::pi_v<traccc::scalar>;
    static constexpr traccc::scalar TWOPI = 2.f * PI;

    // Bring the value within bounds.
    while (phi > PI) {
        phi -= TWOPI;
    }
    while (phi < -PI) {
        phi += TWOPI;
    }
    return phi;
}

inline traccc::scalar angle_mean(traccc::scalar lhs, traccc::scalar rhs) {

    const traccc::scalar diff = wrap_to_pi(lhs - rhs);
    return wrap_to_pi(rhs + 0.5f * diff);
}

}  // namespace

namespace traccc::details {

bool is_same_angle(scalar lhs, scalar rhs, scalar unc) {

    return (std::abs(wrap_to_pi(lhs - rhs)) <=
            (unc * std::abs(angle_mean(lhs, rhs))));
}

}  // namespace traccc::details

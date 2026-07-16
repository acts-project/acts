/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Library include(s).
#include "traccc/performance/details/is_same_scalar.hpp"

// System include(s).
#include <cmath>

namespace traccc::details {

bool is_same_scalar(scalar lhs, scalar rhs, scalar unc) {

    // The difference of the two values is meant to be smaller than
    // their average times the uncertainty.
    return (std::abs(lhs - rhs) <=
            (unc * ((std::abs(lhs) + std::abs(rhs)) / 2.f)));
}

}  // namespace traccc::details

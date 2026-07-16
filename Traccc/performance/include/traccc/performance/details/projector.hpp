/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2025 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Library include(s).
#include "traccc/performance/details/is_same_object.hpp"

// Project include(s).
#include "traccc/definitions/common.hpp"
#include "traccc/definitions/primitives.hpp"

namespace traccc::details {

/// Factory creating instances of "comparator objects" for a given type
///
/// This level of abstraction is necessary to be able to construct comparator
/// objects that would have extra configuration parameters over the reference
/// object and the comparison uncertainty.
///
/// @tparam TYPE The type for which a comparator object should be generated
///
template <typename TYPE>
struct projector {
    static constexpr bool exists = false;
};

template <detray::concepts::algebra algebra_t>
struct projector<traccc::bound_track_parameters<algebra_t>> {
    static constexpr bool exists = true;

    float operator()(const traccc::bound_track_parameters<algebra_t>& i) {
        return static_cast<float>(i.phi());
    }
};

}  // namespace traccc::details

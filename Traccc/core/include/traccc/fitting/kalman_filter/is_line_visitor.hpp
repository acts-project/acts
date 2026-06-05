/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2025 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

#include <detray/geometry/mask.hpp>
#include <detray/geometry/shapes/line.hpp>
#include <detray/geometry/surface.hpp>
#include <detray/utils/type_registry.hpp>

#include "traccc/definitions/qualifiers.hpp"

namespace traccc::detail {

/// @returns true if the surface has "line" shape
template <typename detector_t>
[[nodiscard]] TRACCC_HOST_DEVICE bool constexpr is_line(
    const detray::geometry::surface<detector_t> sf) {
    using algebra_t = typename detector_t::algebra_type;
    using straw_tube = detray::mask<detray::line<false>, algebra_t>;
    using wire_cell = detray::mask<detray::line<true>, algebra_t>;

    using mask_registry_t = typename detector_t::masks;
    if constexpr (detray::types::contains<mask_registry_t, straw_tube> ||
                  detray::types::contains<mask_registry_t, wire_cell>) {
        return (sf.shape_id() ==
                detray::types::id<mask_registry_t, straw_tube>) ||
               (sf.shape_id() == detray::types::id<mask_registry_t, wire_cell>);
    } else {
        return false;
    }
};

}  // namespace traccc::detail

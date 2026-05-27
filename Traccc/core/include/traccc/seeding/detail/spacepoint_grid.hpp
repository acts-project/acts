/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2021-2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "traccc/seeding/grids/axis.hpp"
#include "traccc/seeding/grids/grid2.hpp"
#include "traccc/seeding/grids/populator.hpp"
#include "traccc/seeding/grids/serializer2.hpp"

// Detray include(s)
#include <detray/utils/tuple.hpp>

// VecMem include(s).
#include <vecmem/containers/device_vector.hpp>
#include <vecmem/containers/jagged_device_vector.hpp>
#include <vecmem/containers/jagged_vector.hpp>
#include <vecmem/containers/vector.hpp>

namespace traccc::details {

/// @brief Type definitions for the spacepoint grid.
///
/// This struct contains all the type definitions that are needed to work with
/// the spacepoint grid.
///
struct spacepoint_grid_types {

    /// Spacepoint grid host type
    using host =
        grid2<attach_populator, axis2::circular, axis2::regular, serializer2>;

    /// Spacepoint grid (non-const) device type
    using device =
        grid2<attach_populator, axis2::circular, axis2::regular, serializer2,
              vecmem::device_vector, vecmem::jagged_device_vector>;
    /// Spacepoint grid (const) device type
    using const_device =
        grid2<attach_populator, axis2::circular, axis2::regular, serializer2,
              vecmem::device_vector, vecmem::jagged_device_vector, std::array,
              detray::tuple, const unsigned int>;

    /// Spacepoint grid (non-const) data type
    using data = grid2_data<host>;
    /// Spacepoint grid (const) data type
    using const_data = const_grid2_data<host>;

    /// Spacepoint grid (non-const) view type
    using view = grid2_view<host>;
    /// Spacepoint grid (const) view type
    using const_view = const_grid2_view<host>;

    /// Spacepoint grid buffer type
    using buffer = grid2_buffer<host>;

};  // struct spacepoint_grid_types

}  // namespace traccc::details

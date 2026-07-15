/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "traccc/edm/container.hpp"
#include "traccc/edm/device/doublet_counter.hpp"
#include "traccc/seeding/detail/singlet.hpp"

namespace traccc::device {

/// Doublet of middle-bottom or middle-top spacepoints
struct device_doublet {
    /// bottom (or top) spacepoint location in internal spacepoint container
    sp_location sp2;

    using link_type =
        device::doublet_counter_collection_types::device::size_type;
    /// Link to doublet counter where the middle spacepoint is stored
    link_type counter_link;
};

/// Declare all device doublet collection types
using device_doublet_collection_types = collection_types<device_doublet>;

}  // namespace traccc::device

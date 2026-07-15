/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2021-2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "traccc/seeding/detail/singlet.hpp"

// System include(s).
#include <variant>

namespace traccc {

/// Item: doublet of middle-bottom or middle-top
struct doublet {
    // midle spacepoint location in internal spacepoint container
    sp_location sp1;
    // bottom (or top) spacepoint location in internal spacepoint container
    sp_location sp2;
};

inline TRACCC_HOST_DEVICE bool operator==(const doublet& lhs,
                                          const doublet& rhs) {
    return (lhs.sp1.bin_idx == rhs.sp1.bin_idx &&
            lhs.sp1.sp_idx == rhs.sp1.sp_idx &&
            lhs.sp2.bin_idx == rhs.sp2.bin_idx &&
            lhs.sp2.sp_idx == rhs.sp2.sp_idx);
}

/// Declare all doublet collection types
using doublet_collection_types = collection_types<doublet>;

/// Declare all doublet container types
using doublet_container_types = container_types<std::monostate, doublet>;

}  // namespace traccc

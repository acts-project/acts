/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2021-2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "traccc/definitions/primitives.hpp"
#include "traccc/seeding/detail/doublet.hpp"

// System include(s).
#include <variant>

namespace traccc {

/// Item: triplets of middle-bottom-top
struct triplet {
    // bottom spacepoint location in internal spacepoint container
    sp_location sp1;
    // middle spacepoint location in internal spacepoint container
    sp_location sp2;
    // top spacepoint location in internal spacepoint container
    sp_location sp3;
    // curvtaure of circle estimated from triplet
    scalar curvature;
    // weight of triplet
    scalar weight;
    // z origin of triplet
    scalar z_vertex;
};

inline TRACCC_HOST_DEVICE bool operator==(const triplet& lhs,
                                          const triplet& rhs) {
    return (lhs.sp1.bin_idx == rhs.sp1.bin_idx &&
            lhs.sp1.sp_idx == rhs.sp1.sp_idx &&
            lhs.sp2.bin_idx == rhs.sp2.bin_idx &&
            lhs.sp2.sp_idx == rhs.sp2.sp_idx &&
            lhs.sp3.bin_idx == rhs.sp3.bin_idx &&
            lhs.sp3.sp_idx == rhs.sp3.sp_idx);
}

inline TRACCC_HOST_DEVICE bool operator!=(const triplet& lhs,
                                          const triplet& rhs) {
    return !(lhs == rhs);
}

inline TRACCC_HOST_DEVICE bool operator<(const triplet& lhs,
                                         const triplet& rhs) {
    return lhs.weight < rhs.weight;
}

/// Declare all triplet collection types
using triplet_collection_types = collection_types<triplet>;

/// Declare all triplet container types
using triplet_container_types = container_types<std::monostate, triplet>;

}  // namespace traccc

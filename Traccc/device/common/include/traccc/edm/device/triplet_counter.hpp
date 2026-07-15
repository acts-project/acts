/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2021-2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "traccc/edm/container.hpp"
#include "traccc/seeding/detail/singlet.hpp"

namespace traccc::device {

/// Number of triplets for one specific middle spacepoint.
struct triplet_counter_spM {

    /// Middle spacepoint location in internal spacepoint container
    sp_location spM;

    /// The number of triplets for this middle spacepoint
    unsigned int m_nTriplets = 0;

    /// The position in which these triplets will be added
    unsigned int posTriplets = 0;

};  // struct triplet_counter_spM

/// Declare all triplet counter spM collection types
using triplet_counter_spM_collection_types =
    collection_types<triplet_counter_spM>;

/// Number of triplets for one specific Mid-Bottom Doublet.
struct triplet_counter {

    /// Bottom spacepoint location in internal spacepoint container
    sp_location spB;

    using link_type = triplet_counter_spM_collection_types::device::size_type;
    /// Link to the triplet counter per middle spacepoint
    link_type spM_counter_link;

    /// The number of compatible triplets for this midbot doublet
    unsigned int m_nTriplets = 0;

    /// The position in which these triplets will be added
    unsigned int posTriplets = 0;

};  // struct triplet_counter

/// Declare all triplet counter collection types
using triplet_counter_collection_types = collection_types<triplet_counter>;

}  // namespace traccc::device

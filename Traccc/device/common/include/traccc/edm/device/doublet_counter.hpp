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

/// Number of doublets for one specific middle spacepoint.
struct doublet_counter {

    /// Index of the middle spacepoint.
    sp_location m_spM;

    /// The number of compatible middle-bottom doublets
    unsigned int m_nMidBot = 0;

    /// The number of compatible middle-top doublets
    unsigned int m_nMidTop = 0;

    /// The position in which these middle-bottom doublets will be added
    unsigned int m_posMidBot = 0;

    /// The position in which these middle-top doublets will be added
    unsigned int m_posMidTop = 0;

};  // struct doublet_counter

/// Declare all doublet counter collection types
using doublet_counter_collection_types = collection_types<doublet_counter>;

}  // namespace traccc::device

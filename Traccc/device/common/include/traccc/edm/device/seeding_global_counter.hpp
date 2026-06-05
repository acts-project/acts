/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

namespace traccc::device {

/// Total number of doublets and triplets found in seeding
struct seeding_global_counter {

    /// The total number of middle-bottom doublets
    unsigned int m_nMidBot;

    /// The total number of middle-top doublets
    unsigned int m_nMidTop;

    /// The total number of triplets
    unsigned int m_nTriplets;

};  // struct seeding_global_counter

}  // namespace traccc::device

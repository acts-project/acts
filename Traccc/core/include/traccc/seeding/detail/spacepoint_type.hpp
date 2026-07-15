/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2021-2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

namespace traccc::details {

/// Type of a space point
enum class spacepoint_type : int {
    bottom = 0,  //< The referenced type is a "bottom" spacepoint
    middle = 1,  //< The referenced type is a "middle" spacepoint
    top = 2      //< The referenced type is a "top" spacepoint
};

}  // namespace traccc::details

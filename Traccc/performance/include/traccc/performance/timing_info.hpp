/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// System include(s).
#include <chrono>
#include <iosfwd>
#include <string_view>
#include <utility>
#include <vector>

namespace traccc::performance {

/// Helper type used for timing information storage
using timing_info_pair = std::pair<std::string, std::chrono::nanoseconds>;

/// Struct for storing time measurements collected in timer class
///
struct timing_info {

    /// The low level data.
    std::vector<timing_info_pair> data;

    /// Get the time taken by a given component
    ///
    /// @param timer_name The name of the component
    /// @return The time taken by the component in question
    ///
    std::chrono::nanoseconds get_time(std::string_view timer_name) const;

};  // struct timing_info

/// Printout helper for @c traccc::performance::timing_info
std::ostream& operator<<(std::ostream& out, const timing_info& info);

}  // namespace traccc::performance

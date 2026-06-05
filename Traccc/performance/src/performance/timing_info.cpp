/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Library include(s).
#include "traccc/performance/timing_info.hpp"

// System include(s).
#include <algorithm>
#include <iomanip>
#include <iostream>
#include <stdexcept>

namespace traccc::performance {

std::chrono::nanoseconds timing_info::get_time(
    std::string_view timer_name) const {

    auto it = std::find_if(data.begin(), data.end(),
                           [&timer_name](timing_info_pair itr) {
                               return itr.first == timer_name;
                           });
    if (it == data.end()) {
        throw std::invalid_argument("Unknown component name received");
    }
    return it->second;
}

std::ostream& operator<<(std::ostream& out, const timing_info& info) {

    for (std::size_t i = 0; i < info.data.size(); ++i) {
        const timing_info_pair ti = info.data.at(i);
        out << std::setw(30) << std::right << ti.first << "  "
            << std::chrono::duration_cast<std::chrono::milliseconds>(ti.second)
                   .count()
            << " ms";
        if ((i + 1) < info.data.size()) {
            out << "\n";
        }
    }
    return out;
}

}  // namespace traccc::performance

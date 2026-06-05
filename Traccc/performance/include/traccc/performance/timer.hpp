/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "traccc/performance/timing_info.hpp"

// System include(s).
#include <chrono>
#include <string>
#include <string_view>

namespace traccc::performance {

/// Object used for measuring execution time.
/// Start time measured at construction. End time at destruction.
/// Creating more than one start & stop of timer with the same timer_name &
/// sharing timing info will lead to incrementing the total time for that name
///
class timer {

    public:
    /// Start time measurement
    /// @param timer_name name to be printed out identifying what is measured
    /// @param t_info shared_ptr to timing_info where to store timings
    timer(const std::string_view timer_name, timing_info& t_info);

    /// End time measurement
    ~timer();

    private:
    /// Start time (measured at construct time)
    std::chrono::high_resolution_clock::time_point m_start;

    /// Name of measurement
    std::string m_name;

    /// Shared ptr to timing info where to store elapsed time
    timing_info& m_timing_info;
};  // class timer

}  // namespace traccc::performance

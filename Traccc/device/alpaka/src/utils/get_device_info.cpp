/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2025 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Local include(s).
#include "utils.hpp"

// Project include(s).
#include "traccc/alpaka/utils/get_device_info.hpp"

namespace traccc::alpaka {

std::string get_device_info() {
    int device = 0;
    auto devAcc = ::alpaka::getDevByIdx(::alpaka::Platform<Acc>{}, 0u);
    return std::string("Using Alpaka device: " + ::alpaka::getName(devAcc) +
                       " [id: " + std::to_string(device) + "] ");
}

}  // namespace traccc::alpaka

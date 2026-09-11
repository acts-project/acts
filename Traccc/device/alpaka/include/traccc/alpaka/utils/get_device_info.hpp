/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

#include <string>

namespace traccc::alpaka {

/// Function that prints the current device information to the console.
/// Included as part of the traccc::alpaka namespace, to avoid having to include
/// alpaka headers in any users of the library.
std::string get_device_info();

}  // namespace traccc::alpaka

/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2025 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

namespace traccc {
template <typename T1, typename T2>
struct pair {
    public:
    using first_type = T1;
    using second_type = T2;

    T1 first;
    T2 second;
};
}  // namespace traccc

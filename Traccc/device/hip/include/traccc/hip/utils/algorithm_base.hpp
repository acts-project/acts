/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Local include(s).
#include "traccc/hip/utils/stream.hpp"

// System include(s).
#include <functional>

namespace traccc::hip {

/// Base class for all HIP algorithms
///
/// Holding on to data that all HIP algorithms make use of.
///
class algorithm_base {

    public:
    /// Constructor for the algorithm base
    ///
    /// @param str The HIP stream to perform all operations on
    ///
    explicit algorithm_base(hip::stream& str);

    /// Get the HIP stream of the algorithm
    hip::stream& stream() const;
    /// Get the warp size of the GPU being used
    unsigned int warp_size() const;

    private:
    /// The HIP stream to use
    std::reference_wrapper<hip::stream> m_stream;
    /// Warp size of the GPU being used
    unsigned int m_warp_size;

};  // class algorithm_base

}  // namespace traccc::hip

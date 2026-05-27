/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Local include(s).
#include "traccc/cuda/utils/stream.hpp"

// System include(s).
#include <functional>

namespace traccc::cuda {

/// Base class for all CUDA algorithms
///
/// Holding on to data that all CUDA algorithms make use of.
///
class algorithm_base {

    public:
    /// Constructor for the algorithm base
    ///
    /// @param str The CUDA stream to perform all operations on
    ///
    explicit algorithm_base(cuda::stream& str);

    /// Get the CUDA stream of the algorithm
    cuda::stream& stream() const;
    /// Get the warp size of the GPU being used
    unsigned int warp_size() const;

    private:
    /// The CUDA stream to use
    std::reference_wrapper<cuda::stream> m_stream;
    /// Warp size of the GPU being used
    unsigned int m_warp_size;

};  // class algorithm_base

}  // namespace traccc::cuda

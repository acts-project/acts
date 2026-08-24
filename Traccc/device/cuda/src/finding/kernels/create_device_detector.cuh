/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// CUDA include(s).
#include <cuda_runtime.h>

namespace traccc::cuda {

/// Construct the device detector in place, in global memory
///
/// The detector is built once rather than in every thread, so that the
/// propagation kernels only carry a pointer to it.
///
/// @param stream The stream to launch the kernel on
/// @param in     View of the detector to build from
/// @param out    Global memory to construct the detector into
///
template <typename detector_t>
void create_device_detector(const cudaStream_t& stream,
                            typename detector_t::const_view_type in,
                            detector_t* out);

}  // namespace traccc::cuda

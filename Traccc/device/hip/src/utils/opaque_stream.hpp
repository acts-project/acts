/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2025 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// HIP include(s).
#include <hip/hip_runtime_api.h>

namespace traccc::hip::details {

/// RAII wrapper around @c hipStream_t
///
/// It is used only internally by the HIP library, so it does not need to
/// provide any nice interface.
///
struct opaque_stream {

    /// Default constructor
    opaque_stream(int device);
    /// Destructor
    ~opaque_stream();

    /// Device that the stream is associated to
    int m_device;
    /// Stream managed by the object
    hipStream_t m_stream;

};  // class opaque_stream

}  // namespace traccc::hip::details

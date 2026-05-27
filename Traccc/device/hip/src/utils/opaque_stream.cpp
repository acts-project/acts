/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2025 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Local include(s).
#include "opaque_stream.hpp"

#include "hip_error_handling.hpp"

namespace traccc::hip::details {

opaque_stream::opaque_stream(int device) : m_device{device}, m_stream(nullptr) {

    TRACCC_HIP_ERROR_CHECK(hipStreamCreate(&m_stream));
}

opaque_stream::~opaque_stream() {

    // Don't check the return value of the stream destruction. This is because
    // if the holder of this opaque stream is only destroyed during the
    // termination of the application in which it was created, the HIP runtime
    // may have already deleted all streams by the time that this function would
    // try to delete it.
    //
    // This is not the most robust thing ever, but detecting reliably when this
    // destructor is executed as part of the final operations of an application,
    // would be too platform specific and fragile of an operation.
    [[maybe_unused]] auto code = hipStreamDestroy(m_stream);
}

}  // namespace traccc::hip::details

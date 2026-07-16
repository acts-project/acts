/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2025 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Local include(s).
#include "traccc/hip/utils/stream.hpp"

#include "hip_error_handling.hpp"
#include "opaque_stream.hpp"
#include "utils.hpp"

// HIP include(s).
#include <hip/hip_runtime_api.h>

namespace traccc::hip {

stream::stream(int device) {

    // Make sure that the stream is constructed on the correct device.
    details::select_device dev_selector{
        device == INVALID_DEVICE ? details::get_device() : device};

    // Construct the stream.
    m_stream = std::make_unique<details::opaque_stream>(dev_selector.device());
}

stream::stream(stream&& parent) noexcept = default;

stream::~stream() = default;

stream& stream::operator=(stream&& rhs) noexcept = default;

int stream::device() const {

    return m_stream->m_device;
}

void* stream::hipStream() const {

    return m_stream->m_stream;
}

void stream::synchronize() const {

    TRACCC_HIP_ERROR_CHECK(hipStreamSynchronize(m_stream->m_stream));
}

}  // namespace traccc::hip

/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Local include(s).
#include "traccc/cuda/utils/stream.hpp"

#include "cuda_error_handling.hpp"
#include "get_device.hpp"
#include "select_device.hpp"

// CUDA include(s).
#include <cuda_runtime_api.h>

namespace traccc::cuda {

/// RAII wrapper around @c cudaStream_t
///
/// It is used only internally by the CUDA library, so it does not need to
/// provide any nice interface.
///
struct stream::impl {

    /// Default constructor
    impl(int device) : m_device{device} {
        TRACCC_CUDA_ERROR_CHECK(cudaStreamCreate(&m_stream));
    }
    /// Destructor
    ~impl() {
        // Don't check the return value of the stream destruction. This is
        // because if the holder of this opaque stream is only destroyed during
        // the termination of the application in which it was created, the CUDA
        // runtime may have already deleted all streams by the time that this
        // function would try to delete it.
        //
        // This is not the most robust thing ever, but detecting reliably when
        // this destructor is executed as part of the final operations of an
        // application, would be too platform specific and fragile of an
        // operation.
        cudaStreamDestroy(m_stream);
    }

    /// Device that the stream is associated to
    int m_device;
    /// Stream managed by the object
    cudaStream_t m_stream;

};  // class opaque_stream

stream::stream(int device) {

    // Make sure that the stream is constructed on the correct device.
    details::select_device dev_selector{
        device == INVALID_DEVICE ? details::get_device() : device};

    // Construct the stream.
    m_impl = std::make_unique<impl>(dev_selector.device());
}

stream::stream(stream&& parent) = default;

/// The destructor is implemented explicitly to avoid clients of the class
/// having to know how to destruct @c traccc::cuda::stream::impl objects.
stream::~stream() = default;

stream& stream::operator=(stream&& rhs) = default;

int stream::device() const {

    return m_impl->m_device;
}

void* stream::cudaStream() const {

    return m_impl->m_stream;
}

void stream::synchronize() const {

    TRACCC_CUDA_ERROR_CHECK(cudaStreamSynchronize(m_impl->m_stream));
}

}  // namespace traccc::cuda

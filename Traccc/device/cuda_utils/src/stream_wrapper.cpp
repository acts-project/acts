/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Local include(s).
#include "traccc/cuda/utils/stream_wrapper.hpp"

#include "cuda_error_handling.hpp"

// CUDA include(s).
#include <cuda_runtime_api.h>

namespace traccc::cuda {

stream_wrapper::stream_wrapper(void* stream) : m_stream{stream} {}

int stream_wrapper::device() const {

    // The device ID.
    int device = -1;

#ifdef TRACCC_HAVE_CUDA_STREAM_GET_DEVICE
    // Somewhere around CUDA 12.8 this became available.
    TRACCC_CUDA_ERROR_CHECK(
        cudaStreamGetDevice(static_cast<cudaStream_t>(m_stream), &device));

#else
    // If the CUDA version is too old to support cudaStreamGetDevice, we return
    // the current device instead. This is not ideal, but should be good enough.
    TRACCC_CUDA_ERROR_CHECK(cudaGetDevice(&device));
#endif

    // Return the device ID.
    return device;
}

void* stream_wrapper::cudaStream() const {

    return m_stream;
}

void stream_wrapper::synchronize() const {

    TRACCC_CUDA_ERROR_CHECK(
        cudaStreamSynchronize(static_cast<cudaStream_t>(m_stream)));
}

}  // namespace traccc::cuda

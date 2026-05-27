/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Local include(s).
#include "select_device.hpp"

#include "cuda_error_handling.hpp"

// CUDA include(s).
#include <cuda_runtime_api.h>

namespace traccc::cuda::details {

select_device::select_device(int device) {
    /*
     * When the object is constructed, grab the current device number and
     * store it as a member variable. Then set the device to whatever was
     * specified.
     */
    TRACCC_CUDA_ERROR_CHECK(cudaGetDevice(&m_device));
    TRACCC_CUDA_ERROR_CHECK(cudaSetDevice(device));
}

select_device::~select_device() {
    /*
     * On destruction, reset the device number to whatever it was before the
     * object was constructed.
     */
    TRACCC_CUDA_ERROR_CHECK(cudaSetDevice(m_device));
}

int select_device::device() const {

    return m_device;
}

}  // namespace traccc::cuda::details

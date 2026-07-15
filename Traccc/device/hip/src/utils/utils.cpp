/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2025 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Local include(s).
#include "utils.hpp"

#include "hip_error_handling.hpp"

namespace traccc::hip::details {

int get_device() {

    int d = -1;
    [[maybe_unused]] auto code = hipGetDevice(&d);
    return d;
}

unsigned int get_warp_size(int device) {

    int warp_size = 0;
    TRACCC_HIP_ERROR_CHECK(
        hipDeviceGetAttribute(&warp_size, hipDeviceAttributeWarpSize, device));
    return static_cast<unsigned int>(warp_size);
}

hipStream_t get_stream(const stream& stream) {

    return reinterpret_cast<hipStream_t>(stream.hipStream());
}

select_device::select_device(int device) {
    /*
     * When the object is constructed, grab the current device number and
     * store it as a member variable. Then set the device to whatever was
     * specified.
     */
    TRACCC_HIP_ERROR_CHECK(hipGetDevice(&m_device));
    TRACCC_HIP_ERROR_CHECK(hipSetDevice(device));
}

select_device::~select_device() {
    /*
     * On destruction, reset the device number to whatever it was before the
     * object was constructed.
     */
    TRACCC_HIP_ERROR_CHECK(hipSetDevice(m_device));
}

int select_device::device() const {

    return m_device;
}

}  // namespace traccc::hip::details

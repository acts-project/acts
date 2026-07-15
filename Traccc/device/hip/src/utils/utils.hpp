/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2025 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Local include(s).
#include "traccc/hip/utils/stream.hpp"

// HIP include(s).
#include <hip/hip_runtime_api.h>

namespace traccc::hip::details {

/// Get current HIP device number.
///
/// This function wraps the hipGetDevice function in a way that returns the
/// device number rather than use a reference argument to write to.
///
/// Note that calling the function on a machine with no HIP device does not
/// result in an error, the function just returns -1 in that case.
///
int get_device();

/// Get the warp size for a given device.
///
/// @param device The device to query.
///
/// @return The warp size for the device.
///
unsigned int get_warp_size(int device);

/// Get concrete @c hipStream_t object out of our wrapper
hipStream_t get_stream(const stream& str);

/// Class with RAII mechanism for selecting a HIP device.
///
/// This class can be used to select HIP devices in a modern C++ way, with
/// scope safety. When an object of this class is constructed, it will switch
/// the thread-local device selector to the device number specified in the
/// constructor argument. When this object goes out of scope or gets
/// destructed in any other way, it will restore the device that was set
/// before the object was constructed. This allows us to easily write methods
/// with few side-effects.
///
/// @warning The behaviour of this class is not well-defined if you construct
/// more than one in the same scope.
///
class select_device {

    public:
    /// Constructs the object, switching the current HIP device
    /// to the requested number.
    ///
    /// @param device The HIP device number to switch to.
    ///
    select_device(int device);

    /// Deconstructs the object, returning to the device that was
    /// selected before constructing this object.
    ~select_device();

    /// Return the identifier for the device being seleced
    int device() const;

    private:
    /// The old device number, this is what we restore when the
    /// object goes out of scope.
    int m_device;
};

}  // namespace traccc::hip::details

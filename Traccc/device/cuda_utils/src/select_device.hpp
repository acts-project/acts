/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

namespace traccc::cuda::details {

/// Class with RAII mechanism for selecting a CUDA device.
///
/// This class can be used to select CUDA devices in a modern C++ way, with
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
    /// Constructs the object, switching the current CUDA device
    /// to the requested number.
    ///
    /// @param device The CUDA device number to switch to.
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

}  // namespace traccc::cuda::details

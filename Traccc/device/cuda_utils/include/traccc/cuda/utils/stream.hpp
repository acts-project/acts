/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// System include(s).
#include <memory>

namespace traccc::cuda {

/// Owning wrapper class around @c cudaStream_t
///
/// It is necessary for passing around CUDA stream objects in code that should
/// not be directly exposed to the CUDA header(s).
///
class stream {

    public:
    /// Invalid/default device identifier
    static constexpr int INVALID_DEVICE = -1;

    /// Construct a new stream (possibly for a specified device)
    stream(int device = INVALID_DEVICE);

    /// Move constructor
    stream(stream&& parent);

    /// Destructor
    ~stream();

    /// Move assignment
    stream& operator=(stream&& rhs);

    /// Device that the stream is associated to
    int device() const;

    /// Access a typeless pointer to the managed @c cudaStream_t object
    void* cudaStream() const;

    /// Wait for all queued tasks from the stream to complete
    void synchronize() const;

    private:
    /// Implementation struct hiding the CUDA-specific details
    struct impl;
    /// Smart pointer to the implementation
    std::unique_ptr<impl> m_impl;

};  // class stream

}  // namespace traccc::cuda

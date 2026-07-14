/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2025 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// System include(s).
#include <memory>

namespace traccc::hip {

// Forward declaration(s).
namespace details {
struct opaque_stream;
}

/// Owning wrapper class around @c hipStream_t
///
/// It is necessary for passing around HIP stream objects in code that should
/// not be directly exposed to the HIP header(s).
///
class stream {

    public:
    /// Invalid/default device identifier
    static constexpr int INVALID_DEVICE = -1;

    /// Construct a new stream (possibly for a specified device)
    stream(int device = INVALID_DEVICE);

    /// Move constructor
    stream(stream&& parent) noexcept;

    /// Destructor
    ~stream();

    /// Move assignment
    stream& operator=(stream&& rhs) noexcept;

    /// Device that the stream is associated to
    int device() const;

    /// Access a typeless pointer to the managed @c hipStream_t object
    void* hipStream() const;

    /// Wait for all queued tasks from the stream to complete
    void synchronize() const;

    private:
    /// Smart pointer to the managed @c hipStream_t object
    std::unique_ptr<details::opaque_stream> m_stream;

};  // class stream

}  // namespace traccc::hip

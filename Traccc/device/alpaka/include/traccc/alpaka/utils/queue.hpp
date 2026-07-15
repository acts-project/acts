/**
 * traccc library, part of the ACTS project (R&D line)
 *
 * (c) 2025 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// System include(s).
#include <cstddef>
#include <limits>
#include <memory>

namespace traccc::alpaka {

/// Owning wrapper around @c ::alpaka::Queue
class queue {

    public:
    /// Invalid/default device identifier
    static constexpr std::size_t INVALID_DEVICE =
        std::numeric_limits<std::size_t>::max();

    /// Construct a new stream (possibly for a specified device)
    explicit queue(std::size_t device = INVALID_DEVICE);

    /// Wrap an existing ::alpaka::Queue object
    ///
    /// Without taking ownership of it!
    ///
    explicit queue(void* input_queue);

    /// Move constructor
    queue(queue&& parent) noexcept;

    /// Destructor
    ~queue();

    /// Move assignment
    queue& operator=(queue&& rhs) noexcept;

    /// Access a typeless pointer to the managed @c ::alpaka::Queue object
    void* alpakaQueue();
    /// Access a typeless pointer to the managed @c ::alpaka::Queue object
    const void* alpakaQueue() const;

    /// Wait for all queued tasks from the stream to complete
    void synchronize();

    private:
    /// Type holing the implementation
    struct impl;
    /// Smart pointer to the implementation
    std::unique_ptr<impl> m_impl;

};  // class queue

}  // namespace traccc::alpaka

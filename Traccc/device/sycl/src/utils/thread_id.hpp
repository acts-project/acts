/**
 * traccc library, part of the ACTS project (R&D line)
 *
 * (c) 2024-2025 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "traccc/device/concepts/thread_id.hpp"

// SYCL include(s).
#include <sycl/sycl.hpp>

namespace traccc::sycl::details {

/// A SYCL thread identifier type
template <std::size_t DIMENSIONS>
    requires(DIMENSIONS >= 1 && DIMENSIONS <= 3)
struct thread_id {

    /// Dimensions of the thread identifier
    static constexpr std::size_t dimensions = DIMENSIONS;

    /// Constructor
    ///
    /// @param item The @c ::sycl::nd_item<1> to use for the thread ID.
    ///
    explicit thread_id(const ::sycl::nd_item<dimensions>& item)
        : m_item(item) {}

    /// @name Function(s) implementing @c traccc::device::concepts::thread_id1
    /// @{

    inline unsigned int getLocalThreadId() const {
        return static_cast<unsigned int>(m_item.get_local_linear_id());
    }

    inline unsigned int getLocalThreadIdX() const {
        return static_cast<unsigned int>(m_item.get_local_id(0));
    }

    inline unsigned int getGlobalThreadId() const {
        return static_cast<unsigned int>(m_item.get_global_linear_id());
    }

    inline unsigned int getGlobalThreadIdX() const {
        return static_cast<unsigned int>(m_item.get_global_id(0));
    }

    inline unsigned int getBlockIdX() const {
        return static_cast<unsigned int>(m_item.get_group(0));
    }

    inline unsigned int getBlockDimX() const {
        return static_cast<unsigned int>(m_item.get_local_range(0));
    }

    inline unsigned int getGridDimX() const {
        return static_cast<unsigned int>(m_item.get_global_range(0));
    }

    /// @}

    private:
    /// Item object coming from the SYCL kernel
    const ::sycl::nd_item<dimensions>& m_item;

};  // struct thread_id

/// Template deduction guide for the thread identifier type.
template <int N>
thread_id(::sycl::nd_item<N>) -> thread_id<N>;

/// Verify that @c traccc::sycl::details::thread_id fulfills the
/// @c traccc::device::concepts::thread_id1 concept.
static_assert(traccc::device::concepts::thread_id1<thread_id<1>>);
static_assert(traccc::device::concepts::thread_id1<thread_id<2>>);
static_assert(traccc::device::concepts::thread_id1<thread_id<3>>);

}  // namespace traccc::sycl::details

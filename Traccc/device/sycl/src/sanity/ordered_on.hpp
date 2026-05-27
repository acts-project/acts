/**
 * traccc library, part of the ACTS project (R&D line)
 *
 * (c) 2024-2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "../utils/get_queue.hpp"
#include "traccc/sycl/utils/queue_wrapper.hpp"

// VecMem include(s).
#include <vecmem/memory/memory_resource.hpp>
#include <vecmem/memory/unique_ptr.hpp>
#include <vecmem/utils/copy.hpp>

// SYCL include
#include <sycl/sycl.hpp>

// System include
#include <concepts>

namespace traccc::sycl {

namespace kernels {

/// Kernel used in implementing @c traccc::sycl::is_ordered_on.
///
/// It has to be implemented as a functor, as oneAPI 2024.2.1 gets confused when
/// trying to set up this type of code as a lambda.
///
template <typename CONTAINER, typename R, typename VIEW>
struct is_ordered_on {

    /// Constructor
    is_ordered_on(R relation, VIEW view, bool* out)
        : m_relation(relation), m_view(view), m_out(out) {}

    /// Execution operator for the kernel
    void operator()(::sycl::nd_item<1> item) const {

        std::size_t tid = item.get_global_linear_id();

        const CONTAINER in(m_view);

        if (tid > 0 && tid < in.size()) {
            if (!m_relation(in.at(static_cast<CONTAINER::size_type>(tid - 1)),
                            in.at(static_cast<CONTAINER::size_type>(tid)))) {
                *m_out = false;
            }
        }
    }

    /// The relation object to use
    R m_relation;
    /// View to the input container
    VIEW m_view;
    /// Output boolean
    bool* m_out;
};

}  // namespace kernels

/**
 * @brief Sanity check that a given container is ordered on a given relation.
 *
 * For a container $v$ to be ordered on a relation $R$, it must be the case that
 * for all indices $i$ and $j$, if $i < j$, then $R(i, j)$.
 *
 * @note This function runs in O(n) time.
 *
 * @note Although functions like `std::sort` requires the relation to be strict
 * weak order, this function is more lax in its requirements. Rather, the
 * relation should be a total preorder, i.e. a non-strict weak order.
 *
 * @note For any strict weak order $R$, `is_ordered_on(sort(R, v))` is true.
 *
 * @tparam CONTAINER The type of the (device) container.
 * @tparam R The type of relation $R$, a callable which returns a bool if the
 * first argument can be immediately before the second type.
 * @tparam VIEW The type of the view for the container.
 * @param relation A relation object of type `R`.
 * @param mr A memory resource used for allocating intermediate memory.
 * @param view The container which to check for ordering.
 * @return true If the container is ordered on `R`.
 * @return false Otherwise.
 */
template <typename CONTAINER, std::semiregular R, typename VIEW>
    requires std::regular_invocable<R,
                                    decltype(std::declval<CONTAINER>().at(0)),
                                    decltype(std::declval<CONTAINER>().at(0))>
bool is_ordered_on(R&& relation, vecmem::memory_resource& mr,
                   const vecmem::copy& copy, ::sycl::queue& queue,
                   const VIEW& view) {

    // This should never be a performance-critical step, so we can keep the
    // block size fixed.
    constexpr int block_size = 512;

    // Grab the number of elements in our container.
    const typename VIEW::size_type n = copy.get_size(view);

    // Exit early for empty containers.
    if (n == 0) {
        return true;
    }

    // Initialize the output boolean.
    vecmem::unique_alloc_ptr<bool> out = vecmem::make_unique_alloc<bool>(mr);
    bool initial_out = true;

    ::sycl::event kernel1_memcpy1 =
        queue.memcpy(out.get(), &initial_out, sizeof(bool));

    ::sycl::nd_range<1> kernel_range{
        ::sycl::range<1>(((n + block_size - 1) / block_size) * block_size),
        ::sycl::range<1>(block_size)};

    ::sycl::event kernel1 = queue.submit([&](::sycl::handler& h) {
        h.depends_on(kernel1_memcpy1);
        h.parallel_for<kernels::is_ordered_on<CONTAINER, R, VIEW>>(
            kernel_range, kernels::is_ordered_on<CONTAINER, R, VIEW>(
                              relation, view, out.get()));
    });

    // Copy the output to host, then return it.
    bool host_out;

    queue.memcpy(&host_out, out.get(), sizeof(bool), {kernel1})
        .wait_and_throw();

    return host_out;
}
}  // namespace traccc::sycl

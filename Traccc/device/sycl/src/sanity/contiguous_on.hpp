/**
 * traccc library, part of the ACTS project (R&D line)
 *
 * (c) 2024-2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// VecMem include(s).
#include <vecmem/memory/device_atomic_ref.hpp>
#include <vecmem/memory/memory_resource.hpp>
#include <vecmem/memory/unique_ptr.hpp>
#include <vecmem/utils/copy.hpp>

// SYCL include
#include <sycl/sycl.hpp>

// System include
#include <concepts>
#include <utility>

namespace traccc::sycl {

namespace kernels {

/// Kernel used in implementing @c traccc::sycl::is_contiguous_on
///
/// It has to be implemented as a functor, as oneAPI 2024.2.1 gets confused when
/// trying to set up this type of code as a lambda.
///
template <typename CONTAINER, typename P, typename VIEW,
          std::equality_comparable S>
struct is_contiguous_on_compress_adjacent {

    /// Constructor
    is_contiguous_on_compress_adjacent(P projection, VIEW view,
                                       vecmem::data::vector_view<S> out_view)
        : m_projection(projection), m_view(view), m_out_view(out_view) {}

    /// Execution operator for the kernel
    void operator()(::sycl::nd_item<1> item) const {

        std::size_t tid = item.get_global_linear_id();

        const CONTAINER in(m_view);
        vecmem::device_vector<S> out(m_out_view);

        if (tid > 0 && tid < in.size()) {
            S v1 =
                m_projection(in.at(static_cast<CONTAINER::size_type>(tid - 1)));
            S v2 = m_projection(in.at(static_cast<CONTAINER::size_type>(tid)));

            if (v1 != v2) {
                out.push_back(v2);
            }
        } else if (tid == 0) {
            out.push_back(
                m_projection(in.at(static_cast<CONTAINER::size_type>(tid))));
        }
    }

    /// The projection object to use
    P m_projection;
    /// View to the input container
    VIEW m_view;
    /// Output array
    vecmem::data::vector_view<S> m_out_view;
};

/// Kernel used in implementing @c traccc::sycl::is_contiguous_on
template <typename T>
class is_contiguous_on_all_unique {};

}  // namespace kernels

/**
 * @brief Sanity check that a given container is contiguous on a given
 *        projection.
 *
 * For a container $v$ to be contiguous on a projection $\pi$, it must be the
 * case that for all indices $i$ and $j$, if $v_i = v_j$, then all indices $k$
 * between $i$ and $j$, $v_i = v_j = v_k$.
 *
 * @note This function runs in O(n^2) time.
 *
 * @tparam CONTAINER The type of the (device) container.
 * @tparam P The type of projection $\pi$, a callable which returns some
 * comparable type.
 * @tparam VIEW The type of the view for the container.
 * @param projection A projection object of type `P`.
 * @param mr A memory resource used for allocating intermediate memory.
 * @param view The container which to check for contiguity.
 * @return true If the container is contiguous on `P`.
 * @return false Otherwise.
 */
template <typename CONTAINER, std::semiregular P, typename VIEW>
    requires std::regular_invocable<P,
                                    decltype(std::declval<CONTAINER>().at(0))>
bool is_contiguous_on(P&& projection, vecmem::memory_resource& mr,
                      const vecmem::copy& copy, ::sycl::queue& queue,
                      const VIEW& view) {

    // This should never be a performance-critical step, so we can keep the
    // block size fixed.
    constexpr int local_size = 512;
    constexpr int local_size_2d = 16;

    // Grab the number of elements in our vector.
    const typename VIEW::size_type n = copy.get_size(view);

    // Exit early for empty containers.
    if (n == 0) {
        return true;
    }

    // Get the output type of the projection.
    using projection_t =
        std::invoke_result_t<P, decltype(std::declval<CONTAINER>().at(0))>;

    // Allocate memory for intermediate values and outputs, then set them up.
    vecmem::data::vector_buffer<projection_t> iout(
        n, mr, vecmem::data::buffer_type::resizable);
    copy.setup(iout)->wait();
    vecmem::unique_alloc_ptr<bool> out = vecmem::make_unique_alloc<bool>(mr);

    bool initial_out = true;

    ::sycl::event kernel2_memcpy_evt = queue.copy(&initial_out, out.get(), 1);

    ::sycl::nd_range<1> compress_adjacent_range{
        ::sycl::range<1>(((n + local_size - 1) / local_size) * local_size),
        ::sycl::range<1>(local_size)};

    // Launch the first kernel, which will squash consecutive equal elements
    // into one element.
    queue
        .submit([&](::sycl::handler& h) {
            h.parallel_for<kernels::is_contiguous_on_compress_adjacent<
                CONTAINER, P, VIEW, projection_t>>(
                compress_adjacent_range,
                kernels::is_contiguous_on_compress_adjacent<CONTAINER, P, VIEW,
                                                            projection_t>(
                    projection, view, iout));
        })
        .wait_and_throw();

    typename vecmem::data::vector_view<projection_t>::size_type host_iout_size =
        copy.get_size(iout);
    uint32_t grid_size_rd =
        (host_iout_size + local_size_2d - 1) / local_size_2d;
    ::sycl::nd_range<2> all_unique_range{
        ::sycl::range<2>(grid_size_rd * local_size_2d,
                         grid_size_rd * local_size_2d),
        ::sycl::range<2>(local_size_2d, local_size_2d)};

    // Launch the second kernel, which will check if the values are unique.
    ::sycl::event kernel2_evt = queue.submit([&](::sycl::handler& h) {
        h.depends_on(kernel2_memcpy_evt);
        h.parallel_for<kernels::is_contiguous_on_all_unique<projection_t>>(
            all_unique_range, [in_view = vecmem::get_data(iout),
                               out = out.get()](::sycl::nd_item<2> item) {
                std::size_t tid_x = item.get_global_id(0);
                std::size_t tid_y = item.get_global_id(1);

                const vecmem::device_vector<projection_t> in(in_view);

                if (tid_x < in.size() && tid_y < in.size() && tid_x != tid_y &&
                    in.at(static_cast<CONTAINER::size_type>(tid_x)) ==
                        in.at(static_cast<CONTAINER::size_type>(tid_y))) {
                    *out = false;
                }
            });
    });

    // Get the result from the device and return it.
    bool host_out;

    queue.memcpy(&host_out, out.get(), sizeof(bool), {kernel2_evt})
        .wait_and_throw();

    return host_out;
}
}  // namespace traccc::sycl

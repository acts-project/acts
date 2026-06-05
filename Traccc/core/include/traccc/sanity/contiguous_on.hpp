/**
 * traccc library, part of the ACTS project (R&D line)
 *
 * (c) 2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// System include
#include <concepts>
#include <memory>
#include <unordered_set>
#include <utility>

namespace traccc::host {

/**
 * @brief Sanity check that a given container is contiguous on a given
 *        projection.
 *
 * For a container $c$ to be contiguous on a projection $\pi$, it must be the
 * case that for all indices $i$ and $j$, if $v_i = v_j$, then all indices $k$
 * between $i$ and $j$, $v_i = v_j = v_k$.
 *
 * @note This function runs in O(n^2) time.
 *
 * @tparam P The type of projection $\pi$, a callable which returns some
 * comparable type.
 * @tparam CONTAINER The (SoA) container type
 * @param projection A projection object of type `P`.
 * @param mr A memory resource used for allocating intermediate memory.
 * @param vector The vector which to check for contiguity.
 * @return true If the vector is contiguous on `P`.
 * @return false Otherwise.
 */
template <std::semiregular P, typename CONTAINER>
    requires std::regular_invocable<P,
                                    decltype(std::declval<CONTAINER>().at(0))>
bool is_contiguous_on(P&& projection, const CONTAINER& in) {

    // Grab the number of elements in our container.
    typename CONTAINER::size_type n = in.size();

    // Get the output type of the projection.
    using projection_t =
        std::invoke_result_t<P, decltype(std::declval<CONTAINER>().at(0))>;

    // Allocate memory for intermediate values and outputs, then set them up.
    std::unique_ptr<projection_t[]> iout = std::make_unique<projection_t[]>(n);
    std::size_t iout_size = 0;

    // Compress adjacent elements
    for (std::size_t i = 0; i < n; ++i) {
        if (i == 0) {
            iout[iout_size++] = projection(
                in.at(static_cast<typename CONTAINER::size_type>(i)));
        } else {
            projection_t v = projection(
                in.at(static_cast<typename CONTAINER::size_type>(i)));

            if (v != iout[iout_size - 1]) {
                iout[iout_size++] = v;
            }
        }
    }

    // Check whether all elements are unique
    std::unordered_set<projection_t> seen;

    for (std::size_t i = 0; i < iout_size; ++i) {
        projection_t& v = iout[i];

        if (seen.count(v) == 1) {
            return false;
        } else {
            seen.insert(v);
        }
    }

    return true;
}

}  // namespace traccc::host

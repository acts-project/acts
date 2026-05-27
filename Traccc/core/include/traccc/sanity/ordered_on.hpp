/**
 * traccc library, part of the ACTS project (R&D line)
 *
 * (c) 2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// VecMem include(s).
#include <vecmem/containers/data/vector_view.hpp>
#include <vecmem/containers/device_vector.hpp>
#include <vecmem/memory/memory_resource.hpp>
#include <vecmem/utils/copy.hpp>

// System include
#include <concepts>
#include <utility>

namespace traccc::host {
/**
 * @brief Sanity check that a given vector is ordered on a given relation.
 *
 * For a vector $v$ to be ordered on a relation $R$, it must be the case that
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
 * @tparam R The type of relation $R$, a callable which returns a bool if the
 * first argument can be immediately before the second type.
 * @tparam T The type of the vector.
 * @param relation A relation object of type `R`.
 * @param mr A memory resource used for allocating intermediate memory.
 * @param vector The vector which to check for ordering.
 * @return true If the vector is ordered on `R`.
 * @return false Otherwise.
 */
template <std::semiregular R, typename CONTAINER>
    requires std::regular_invocable<R,
                                    decltype(std::declval<CONTAINER>().at(0)),
                                    decltype(std::declval<CONTAINER>().at(0))>
bool is_ordered_on(R&& relation, const CONTAINER& in) {

    // Grab the number of elements in our vector.
    typename CONTAINER::size_type n = in.size();

    // Check for orderedness.
    for (typename CONTAINER::size_type i = 1; i < n; ++i) {
        if (!relation(in.at(i - 1), in.at(i))) {
            return false;
        }
    }

    return true;
}
}  // namespace traccc::host

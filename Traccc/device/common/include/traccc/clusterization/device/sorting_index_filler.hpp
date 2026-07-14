/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2025 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "traccc/definitions/qualifiers.hpp"

// VecMem include(s).
#include <vecmem/containers/data/vector_view.hpp>

namespace traccc::device {

/// Helper functor for filling the indices used during measurement sorting.
class sorting_index_filler {

    public:
    /// Constructor with the indices view to fill.
    ///
    /// @param indices_view The view of the indices to fill
    ///
    explicit TRACCC_HOST_DEVICE sorting_index_filler(
        vecmem::data::vector_view<unsigned int> indices)
        : m_indices{indices} {}

    /// The operator filling the indices.
    ///
    /// @param index The index to fill
    ///
    TRACCC_HOST_DEVICE void operator()(unsigned int& index) const {
        index = static_cast<unsigned int>(&index - m_indices.ptr());
    }

    private:
    /// The view of the indices to fill.
    vecmem::data::vector_view<unsigned int> m_indices;

};  // class sorting_index_filler

}  // namespace traccc::device

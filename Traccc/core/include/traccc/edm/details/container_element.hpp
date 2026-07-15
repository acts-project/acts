/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2021-2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "traccc/definitions/qualifiers.hpp"

namespace traccc {

/**
 * @brief View class for an element in a header-vector container.
 *
 * In order to enforce certain invariants on the @c container header-vector
 * type, we access elements (which, of course, are product types of a header
 * and a vector) through this wrapper class. This class provides
 * low-overhead access to the struct-of-arrays container class to emulate an
 * array-of-structs architecture.
 *
 * @tparam header_reference_t The type of reference of the header object.
 * @tparam vector_reference_t The fully qualified vector reference type.
 */
template <typename header_reference_t, typename vector_reference_t>
struct container_element {

    /// Header reference type
    using header_reference = header_reference_t;
    /// Vector reference type
    using vector_reference = vector_reference_t;

    /**
     * @brief Construct a new container element view.
     *
     * This constructor is extremely trivial, as it simply takes a reference
     * to a header, a reference to a vector, and saves them in this object's
     * internal state.
     *
     * @param[in] h The header object reference.
     * @param[in] v The vector object reference.
     */
    TRACCC_HOST_DEVICE
    container_element(header_reference h, vector_reference v)
        : header(h), items(v) {}

    header_reference header;
    vector_reference items;
};  // struct container_element

}  // namespace traccc

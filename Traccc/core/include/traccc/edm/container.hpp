/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2021-2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "traccc/definitions/qualifiers.hpp"
#include "traccc/edm/details/device_container.hpp"
#include "traccc/edm/details/host_container.hpp"
#include "traccc/utils/type_traits.hpp"

// VecMem include(s).
#include <vecmem/containers/data/jagged_vector_buffer.hpp>
#include <vecmem/containers/data/jagged_vector_data.hpp>
#include <vecmem/containers/data/jagged_vector_view.hpp>
#include <vecmem/containers/data/vector_buffer.hpp>
#include <vecmem/containers/data/vector_view.hpp>
#include <vecmem/containers/device_vector.hpp>
#include <vecmem/containers/vector.hpp>

// System include(s).
#include <type_traits>

namespace traccc {

/// @name Types used to send data back and forth between host and device code
/// @{

/// Structure holding (some of the) data about the container in host code
template <typename header_t, typename item_t>
struct container_data {
    using header_vector = vecmem::data::vector_view<header_t>;
    using item_vector = vecmem::data::jagged_vector_data<item_t>;
    header_vector headers;
    item_vector items;
};

/// Structure holding (all of the) data about the container in host code
template <typename header_t, typename item_t>
struct container_buffer {
    using header_vector = vecmem::data::vector_buffer<header_t>;
    using item_vector = vecmem::data::jagged_vector_buffer<item_t>;
    header_vector headers;
    item_vector items;
};

/// Structure used to send the data about the container to device code
///
/// This is the type that can be passed to device code as-is. But since in
/// host code one needs to manage the data describing a
/// @c traccc::container either using @c traccc::container_data or
/// @c traccc::container_buffer, it needs to have constructors from
/// both of those types.
///
/// In fact it needs to be created from one of those types, as such an
/// object can only function if an instance of one of those types exists
/// alongside it as well.
///
template <typename header_t, typename item_t>
struct container_view {

    /// Type for the header vector (view)
    using header_vector = vecmem::data::vector_view<header_t>;
    /// Type for the item vector (view)
    using item_vector = vecmem::data::jagged_vector_view<item_t>;

    /// Constructor from a @c container_data object
    template <
        typename other_header_t, typename other_item_t,
        std::enable_if_t<details::is_same_nc<header_t, other_header_t>::value,
                         bool> = true,
        std::enable_if_t<details::is_same_nc<item_t, other_item_t>::value,
                         bool> = true>
    container_view(const container_data<other_header_t, other_item_t>& data)
        : headers(data.headers), items(data.items) {}

    /// Constructor from a @c container_buffer object
    template <
        typename other_header_t, typename other_item_t,
        std::enable_if_t<details::is_same_nc<header_t, other_header_t>::value,
                         bool> = true,
        std::enable_if_t<details::is_same_nc<item_t, other_item_t>::value,
                         bool> = true>
    container_view(const container_buffer<other_header_t, other_item_t>& buffer)
        : headers(buffer.headers), items(buffer.items) {}

    /// Constructor from a non-const view
    template <
        typename other_header_t, typename other_item_t,
        std::enable_if_t<details::is_same_nc<header_t, other_header_t>::value,
                         bool> = true,
        std::enable_if_t<details::is_same_nc<item_t, other_item_t>::value,
                         bool> = true>
    container_view(const container_view<other_header_t, other_item_t>& parent)
        : headers(parent.headers), items(parent.items) {}

    /// View of the data describing the headers
    header_vector headers;

    /// View of the data describing the items
    item_vector items;
};

/// Helper function for making a "simple" object out of the container
/// (non-const)
template <typename header_t, typename item_t>
inline container_data<header_t, item_t> get_data(
    host_container<header_t, item_t>& cc,
    vecmem::memory_resource* resource = nullptr) {
    return {{vecmem::get_data(cc.get_headers())},
            {vecmem::get_data(cc.get_items(), resource)}};
}

/// Helper function for making a "simple" object out of the container (const)
template <typename header_t, typename item_t>
inline container_data<const header_t, const item_t> get_data(
    const host_container<header_t, item_t>& cc,
    vecmem::memory_resource* resource = nullptr) {
    return {{vecmem::get_data(cc.get_headers())},
            {vecmem::get_data(cc.get_items(), resource)}};
}

/// Type trait defining all "collection types" for an EDM class
template <typename item_t>
struct collection_types {

    /// @c item_t must not be a constant type
    static_assert(std::is_const<item_t>::value == false,
                  "The template parameter must not be a constant type");

    /// Host collection for @c item_t
    using host = vecmem::vector<item_t>;
    /// Non-const device collection for @c item_t
    using device = vecmem::device_vector<item_t>;
    /// Constant device collection for @c item_t
    using const_device = vecmem::device_vector<const item_t>;

    /// Non-constant view of an @c item_t collection
    using view = vecmem::data::vector_view<item_t>;
    /// Constant view of an @c item_t collection
    using const_view = vecmem::data::vector_view<const item_t>;

    /// Buffer for an @c item_t collection
    using buffer = vecmem::data::vector_buffer<item_t>;

};  // struct collection_types

/// Type trait defining all "container types" for an EDM class pair
template <typename header_t, typename item_t>
struct container_types {

    /// @c header_t must not be a constant type
    static_assert(std::is_const<header_t>::value == false,
                  "The header type must not be constant");
    /// @c item_t must not be a constant type
    static_assert(std::is_const<item_t>::value == false,
                  "The item type must not be constant");

    /// Host container for @c header_t and @c item_t
    using host = host_container<header_t, item_t>;
    /// Non-const device container for @c header_t and @c item_t
    using device = device_container<header_t, item_t>;
    /// Constant device container for @c header_t and @c item_t
    using const_device = device_container<const header_t, const item_t>;

    /// Non-constant view of an @c header_t / @c item_t container
    using view = container_view<header_t, item_t>;
    /// Constant view of an @c header_t / @c item_t container
    using const_view = container_view<const header_t, const item_t>;

    /// Non-constant data for an @c header_t / @c item_t container
    using data = container_data<header_t, item_t>;
    /// Constant data for an @c header_t / @c item_t container
    using const_data = container_data<const header_t, const item_t>;

    /// Buffer for an @c header_t / @c item_t container
    using buffer = container_buffer<header_t, item_t>;

};  // struct container_types

}  // namespace traccc

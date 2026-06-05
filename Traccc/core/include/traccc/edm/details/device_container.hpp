/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2021-2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "traccc/definitions/qualifiers.hpp"
#include "traccc/edm/details/container_base.hpp"

// VecMem include(s).
#include <vecmem/containers/device_vector.hpp>
#include <vecmem/containers/jagged_device_vector.hpp>

// System include(s).
#include <cassert>

namespace traccc {

/// Host container describing objects in a given event
template <typename header_t, typename item_t>
class device_container
    : public container_base<header_t, item_t, vecmem::device_vector,
                            vecmem::jagged_device_vector> {
    public:
    /// Base class type
    using base_type = container_base<header_t, item_t, vecmem::device_vector,
                                     vecmem::jagged_device_vector>;

    /// Inherit all of the base class's constructors
    using base_type::base_type;

    /// Constructor from a "view object"
    ///
    /// Note that it may also be a "data" or "buffer" object. It just needs to
    /// look like a "view object".
    ///
    template <template <typename, typename> class view_t>
    TRACCC_HOST_DEVICE device_container(const view_t<header_t, item_t>& view)
        : base_type(view.headers, view.items) {}

    /**
     * @brief Bounds-checking mutable element accessor.
     */
    TRACCC_HOST_DEVICE
    typename base_type::element_view at(typename base_type::size_type i) {
        return {this->m_headers.at(i), this->m_items.at(i)};
    }

    /**
     * @brief Bounds-checking immutable element accessor.
     */
    TRACCC_HOST_DEVICE
    typename base_type::const_element_view at(
        typename base_type::size_type i) const {
        return {this->m_headers.at(i), this->m_items.at(i)};
    }

    /**
     * @brief Bounds-checking mutable item vector element accessor.
     */
    TRACCC_HOST_DEVICE
    typename base_type::item_type& at(
        const typename base_type::link_type& link) {
        return this->m_items.at(link.first).at(link.second);
    }

    /**
     * @brief Bounds-checking immutable item vector element accessor.
     */
    TRACCC_HOST_DEVICE
    const typename base_type::item_type& at(
        const typename base_type::link_type& link) const {
        return this->m_items.at(link.first).at(link.second);
    }

    /**
     * @brief Mutable element accessor.
     */
    TRACCC_HOST_DEVICE
    typename base_type::element_view operator[](
        typename base_type::size_type i) {
        return {this->m_headers[i], this->m_items[i]};
    }

    /**
     * @brief Immutable element accessor.
     */
    TRACCC_HOST_DEVICE
    typename base_type::const_element_view operator[](
        typename base_type::size_type i) const {
        return {this->m_headers[i], this->m_items[i]};
    }

    /**
     * @brief Mutable item vector element accessor.
     */
    TRACCC_HOST_DEVICE
    typename base_type::item_type& operator[](
        const typename base_type::link_type& link) {
        return this->m_items[link.first][link.second];
    }

    /**
     * @brief Immutable item vector element accessor.
     */
    TRACCC_HOST_DEVICE
    const typename base_type::item_type& operator[](
        const typename base_type::link_type& link) const {
        return this->m_items[link.first][link.second];
    }

    /**
     * @brief Return the size of the container.
     *
     * In principle, the size of the two internal vectors should always be
     * equal, but we can assert this at runtime for debug builds.
     */
    TRACCC_HOST_DEVICE
    typename base_type::size_type size(void) const {
        assert(this->m_headers.size() == this->m_items.size());
        return this->m_headers.size();
    }

    /**
     * @breif Get number of items of jagged vector
     */
    TRACCC_HOST_DEVICE
    typename base_type::size_type total_size() const {
        typename base_type::size_type ret = 0;
        for (const typename base_type::item_vector::value_type& item :
             this->m_items) {
            ret += item.size();
        }
        return ret;
    }

};  // class device_container

}  // namespace traccc

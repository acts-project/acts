/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2021-2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "traccc/definitions/qualifiers.hpp"
#include "traccc/edm/details/container_element.hpp"
#include "traccc/utils/pair.hpp"

// VecMem include(s).
#include <vecmem/memory/memory_resource.hpp>

// System include(s).
#include <cassert>
#include <stdexcept>
#include <type_traits>

namespace traccc {

/// Container base class for describing objects in a given event
///
/// This is the generic container of the code, holding all relevant
/// information about objcts in a given event.
///
/// It can be instantiated with different vector types, to be able to use
/// the same container type in both host and device code.
///
/// It also can be instantiated with different edm types represented by
/// header and item type.
///
template <typename header_t, typename item_t,
          template <typename> class vector_t,
          template <typename> class jagged_vector_t,
          template <typename, typename> class pair_t = traccc::pair>
class container_base {
    public:
    /// @name Type definitions
    /// @{

    /// Header type
    using header_type = header_t;

    /// Item type
    using item_type = item_t;

    /// Vector type used by the container
    template <typename T>
    using vector_type = vector_t<T>;

    /// Jagged vector type used by the container
    template <typename T>
    using jagged_vector_type = jagged_vector_t<T>;

    /// The header vector type
    using header_vector = vector_type<header_type>;

    /// The item vector type
    using item_vector = jagged_vector_type<item_type>;

    /**
     * @brief The size type of this container, which is the type by which
     * its elements are indexed.
     */
    using size_type = typename header_vector::size_type;

    /// The element link type
    using link_type = pair_t<typename header_vector::size_type,
                             typename item_vector::size_type>;

    /// @}

    /**
     * @brief The type name of the element view which is returned by various
     * methods in this class.
     */
    using element_view = container_element<typename header_vector::reference,
                                           typename item_vector::reference>;

    /**
     * @brief The type name of the constant element view which is returned
     * by various methods in this class.
     */
    using const_element_view =
        container_element<typename header_vector::const_reference,
                          typename item_vector::const_reference>;

    /**
     * We need to assert that the header vector and the outer layer of the
     * jagged vector have the same size type, so they can be indexed using
     * the same type.
     */
    static_assert(
        std::is_convertible<typename header_vector::size_type,
                            typename item_vector::size_type>::value,
        "Size type for container header and item vectors must be the same.");

    /// Default copy-constructor
    container_base(const container_base&) = default;
    /// Default move-constructor
    container_base(container_base&&) = default;

    /**
     * @brief (Copy) Constructor from a header and item vector
     *
     * To enforce the invariant that both vectors must be the same size, we
     * check this in the constructor. This is also checked in release
     * builds.
     */
    TRACCC_HOST
    container_base(const header_vector& hv, const item_vector& iv)
        : m_headers(hv), m_items(iv) {
        if (m_headers.size() != m_items.size()) {
            throw std::logic_error("Header and item length not equal.");
        }
    }

    /**
     * @brief (Move) Constructor from a header and item vector
     *
     * To enforce the invariant that both vectors must be the same size, we
     * check this in the constructor. This is also checked in release
     * builds.
     */
    TRACCC_HOST
    container_base(header_vector&& hv, item_vector&& iv)
        : m_headers(hv), m_items(iv) {
        if (m_headers.size() != m_items.size()) {
            throw std::logic_error("Header and item length not equal.");
        }
    }

    /**
     * @brief Constructor from a pair of "view type" objects
     */
    template <typename header_vector_tp, typename item_vector_tp>
        requires(std::constructible_from<header_vector,
                                         const header_vector_tp&> &&
                 std::constructible_from<item_vector, const item_vector_tp&>)
    TRACCC_HOST_DEVICE container_base(const header_vector_tp& hv,
                                      const item_vector_tp& iv)
        : m_headers(hv), m_items(iv) {

        assert(m_headers.size() == m_items.size());
    }

    /**
     * @brief Constructor with vector size and memory resource .
     */
    template <typename size_type>
    TRACCC_HOST explicit container_base(size_type size,
                                        vecmem::memory_resource* mr)
        : m_headers(size, mr), m_items(size, mr) {}

    /**
     * @brief Constructor with memory resource .
     */

    TRACCC_HOST explicit container_base(vecmem::memory_resource* mr)
        : m_headers(mr), m_items(mr) {}

    /**
     * @brief Default Constructor
     */
    container_base() = default;

    /// Default (copy) assignment operator
    container_base& operator=(const container_base&) = default;
    /// Default (move) assignment operator
    container_base& operator=(container_base&&) = default;

    /**
     * @brief Accessor method for the internal header vector.
     */
    TRACCC_HOST_DEVICE
    const header_vector& get_headers() const { return m_headers; }

    /**
     * @brief Non-const accessor method for the internal header vector.
     *
     * @warning Do not use this function! It is dangerous, and risks breaking
     * invariants!
     */
    TRACCC_HOST_DEVICE
    header_vector& get_headers() { return m_headers; }

    /**
     * @brief Accessor method for the internal item vector-of-vectors.
     */
    TRACCC_HOST_DEVICE
    const item_vector& get_items() const { return m_items; }

    /**
     * @brief Non-const accessor method for the internal item vector-of-vectors.
     *
     * @warning Do not use this function! It is dangerous, and risks breaking
     * invariants!
     */
    TRACCC_HOST_DEVICE
    item_vector& get_items() { return m_items; }

    protected:
    /// Headers information related to the objects in the event
    header_vector m_headers;

    /// All objects in the event
    item_vector m_items;

};  // class container_base

}  // namespace traccc

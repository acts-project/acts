/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2020-2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "traccc/definitions/qualifiers.hpp"

// Detray include(s)
#include <detray/definitions/algorithms.hpp>
#include <detray/utils/invalid_values.hpp>

// VecMem include(s).
#include <vecmem/containers/data/jagged_vector_buffer.hpp>
#include <vecmem/containers/data/jagged_vector_view.hpp>
#include <vecmem/containers/data/vector_buffer.hpp>
#include <vecmem/containers/data/vector_view.hpp>
#include <vecmem/containers/device_vector.hpp>
#include <vecmem/containers/jagged_device_vector.hpp>

// System include(s).
#include <algorithm>
#include <limits>

namespace traccc {

/** A replace populator that swaps whatever current value in the
 * bin with the new one.
 *
 * @tparam value_t the type of a single stored object
 *
 * @note bare_value and store_value are identicial in this case
 **/
template <template <typename...> class vector_t = vecmem::vector,
          template <typename...> class jagged_vector_t = vecmem::jagged_vector,
          template <typename, std::size_t> class array_t = std::array,
          typename value_t = unsigned int, bool kSORT = false,
          unsigned int kDIM = 1u>
struct replace_populator {
    DETRAY_HOST_DEVICE
    explicit replace_populator(
        const value_t invalid = detray::detail::invalid_value<value_t>())
        : m_invalid(invalid) {}

    value_t m_invalid;

    using bare_value = value_t;
    using store_value = value_t;
    using serialized_storage = vector_t<store_value>;

    using vector_view_type = vecmem::data::vector_view<store_value>;
    using const_vector_view_type = vecmem::data::vector_view<const store_value>;
    using vector_data_type = vecmem::data::vector_view<store_value>;
    using const_vector_data_type = vecmem::data::vector_view<const store_value>;
    using vector_buffer_type = vecmem::data::vector_buffer<store_value>;
    using buffer_size_type = typename vector_view_type::size_type;

    /** Swap the stored value with a new bare value
     *
     * @param stored the stored value for the population
     * @param bvalue the new value to be added
     **/
    DETRAY_HOST_DEVICE
    void operator()(store_value &stored, bare_value &&bvalue) const {
        stored = std::move(bvalue);
    }

    /** Create a sequence of bare values, independent of the store_value.
     *
     * @param stored the stored value
     *
     * @return a sequence of bare values
     */
    DETRAY_HOST_DEVICE
    vector_t<bare_value> sequence(store_value &stored) const {
        if (stored != m_invalid) {
            return {stored};
        }
        return {};
    }

    /** Shift operation for unified memory block
     *
     * @param stored the stored value
     * @param offset is the shift offset
     *
     **/
    DETRAY_HOST_DEVICE
    void shift(store_value &stored, const bare_value &offset) const {
        stored += offset;
    }

    /** Return an initialized bin value
     */
    DETRAY_HOST_DEVICE
    store_value init() const { return m_invalid; }
};

/** A complete populator that adds values to the internal
 * store array until it is completed, ignored afterwards.
 *
 * @tparam kDIM the dimension of the underlying stored array
 * @tparam kSORT a sorting flag
 * @tparam value_t the type of a single stored object
 * @tparam m_invalid the chosen invalid type
 *
 * @note bare_value and store_value are different in this case
 **/

template <template <typename...> class vector_t = vecmem::vector,
          template <typename...> class jagged_vector_t = vecmem::jagged_vector,
          template <typename, std::size_t> class array_t = std::array,
          typename value_t = unsigned int, bool kSORT = false,
          unsigned int kDIM = 1u>
struct complete_populator {
    DETRAY_HOST_DEVICE
    explicit complete_populator(
        const value_t invalid = detray::detail::invalid_value<value_t>())
        : m_invalid(invalid) {}

    value_t m_invalid;

    using bare_value = value_t;
    using store_value = array_t<bare_value, kDIM>;
    using serialized_storage = vector_t<store_value>;

    using vector_view_type = vecmem::data::vector_view<store_value>;
    using const_vector_view_type = vecmem::data::vector_view<const store_value>;
    using vector_data_type = vecmem::data::vector_view<store_value>;
    using const_vector_data_type = vecmem::data::vector_view<const store_value>;
    using vector_buffer_type = vecmem::data::vector_buffer<store_value>;
    using buffer_size_type = typename vector_view_type::size_type;

    /** Complete the stored value with a new bare value - for host
     *
     * @param stored the stored value for the population
     * @param bvalue the new value to be added
     **/
    DETRAY_HOST_DEVICE
    void operator()(store_value &stored, bare_value &&bvalue) const {
        for (auto &val : stored) {
            if (val == m_invalid) {
                val = std::move(bvalue);
                break;
            }
        }
        if constexpr (kSORT) {
            detray::sequential_sort(stored.begin(), stored.end());
        }
    }

    /** Create a sequence of bare values, independent of the store_value.
     *
     * @param stored the stored value
     *
     * @return a sequence of bare values, @note it will ignore invalid entries
     */
    DETRAY_HOST_DEVICE
    vector_t<bare_value> sequence(store_value &stored) const {
        vector_t<bare_value> s;
        s.reserve(kDIM);
        for (const auto &val : stored) {
            if (val != m_invalid) {
                s.push_back(val);
            }
        }
        return s;
    }

    /** Shift operation for unified memory block
     *
     * @param stored the stored value
     * @param offset is the shift offset
     *
     **/
    DETRAY_HOST_DEVICE
    void shift(store_value &stored, const bare_value &offset) const {
        std::ranges::for_each(stored, [&](auto &d) { d += offset; });
    }

    /** Return an initialized bin value
     **/
    DETRAY_HOST_DEVICE
    store_value init() const {

        store_value init_bin;
        for (auto &val : init_bin) {
            val = m_invalid;
        }
        return init_bin;
    }
};

/** An attach populator that adds the new value to the
 *
 * @tparam kSORT the sorting directive
 * @tparam value_t the type of a single stored object
 *
 * @note bare_value and store_value are identicial in this case
 **/
template <template <typename...> class vector_t = vecmem::vector,
          template <typename...> class jagged_vector_t = vecmem::jagged_vector,
          template <typename, std::size_t> class array_t = std::array,
          typename value_t = unsigned int, bool kSORT = false,
          unsigned int kDIM = 1u>
struct attach_populator {
    DETRAY_HOST_DEVICE
    explicit attach_populator(
        const value_t invalid = detray::detail::invalid_value<value_t>())
        : m_invalid(invalid) {}

    value_t m_invalid;

    using bare_value = value_t;
    using store_value = vector_t<bare_value>;
    using serialized_storage = jagged_vector_t<bare_value>;

    using vector_view_type = vecmem::data::jagged_vector_view<bare_value>;
    using const_vector_view_type =
        vecmem::data::jagged_vector_view<const bare_value>;
    using vector_data_type = vecmem::data::jagged_vector_data<bare_value>;
    using const_vector_data_type =
        vecmem::data::jagged_vector_data<const bare_value>;
    using vector_buffer_type = vecmem::data::jagged_vector_buffer<bare_value>;
    using buffer_size_type = std::vector<typename vector_view_type::size_type>;

    /** Add a new value to the stored value - for host vector
     *
     * @param stored the stored value for the population
     * @param bvalue the new value to be added
     **/
    DETRAY_HOST
    void operator()(store_value &stored, bare_value &&bvalue) const {
        stored.push_back(std::move(bvalue));
        if (kSORT) {
            std::ranges::sort(stored);
        }
    }

    /** Add a new value to the stored value - for device vector
     *
     * @param stored the stored value for the population
     * @param bvalue the new value to be added
     **/
#if defined(__CUDACC__)  // to resolve ambiguoty from host side
    DETRAY_DEVICE
    void operator()(store_value stored, bare_value &&bvalue) const {
        stored.push_back(std::move(bvalue));
    }
#endif

    /** Create a sequence of bare values, independent of the store_value.
     *
     * @param stored the stored value
     *
     * @return a sequence of bare values
     */
    DETRAY_HOST_DEVICE
    vector_t<bare_value> sequence(store_value &stored) const { return stored; }

    /** Shift operation for unified memory block
     *
     * @param stored the stored value
     * @param offset is the shift offset
     *
     **/
    DETRAY_HOST_DEVICE
    void shift(store_value &stored, const bare_value &offset) const {
        std::ranges::for_each(stored, [&](auto &d) { d += offset; });
    }

    /** Return an initialized bin value
     **/
    DETRAY_HOST_DEVICE
    store_value init() const { return {}; }
};

}  // namespace traccc

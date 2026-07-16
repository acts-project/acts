/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2020-2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "traccc/definitions/primitives.hpp"
#include "traccc/definitions/qualifiers.hpp"
#include "traccc/seeding/grids/axis.hpp"

// Detray include(s).
#include <detray/utils/concepts.hpp>
#include <detray/utils/invalid_values.hpp>
#include <detray/utils/tuple.hpp>

// VecMem include(s).
#include <vecmem/containers/data/buffer_type.hpp>
#include <vecmem/memory/memory_resource.hpp>

// System include(s)
#include <algorithm>

namespace traccc {

/** A two-dimensional grid for object storage
 *
 * @tparam populator_t  is a prescription what to do when a bin gets
 * populated, it broadcasts also the value type
 * @tparam tparam axis_p0_t the type of the first axis
 * @tparam tparam axis_p1_t the type of the second axis
 * @tparam serialzier_t  type of the serializer to the storage represenations
 *
 **/
template <template <template <typename...> class, template <typename...> class,
                    template <typename, std::size_t> class, typename, bool,
                    unsigned int>
          class populator_t,
          template <template <typename, std::size_t> class,
                    template <typename...> class>
          class axis_p0_t,
          template <template <typename, std::size_t> class,
                    template <typename...> class>
          class axis_p1_t,
          typename serializer_t,
          template <typename...> class vector_t = vecmem::vector,
          template <typename...> class jagged_vector_t = vecmem::jagged_vector,
          template <typename, std::size_t> class array_t = std::array,
          template <typename...> class tuple_t = detray::tuple,
          typename value_t = unsigned int, bool kSORT = false,
          unsigned int kDIM = 1u>
class grid2 {

    public:
    using populator_type =
        populator_t<vector_t, jagged_vector_t, array_t, value_t, kSORT, kDIM>;
    using serializer_type = serializer_t;
    using axis_p0_type = axis_p0_t<array_t, vector_t>;
    using axis_p1_type = axis_p1_t<array_t, vector_t>;
    using bare_value = typename populator_type::bare_value;
    using serialized_storage = typename populator_type::serialized_storage;

    template <typename neighbor_t>
    using neighborhood = array_t<array_t<neighbor_t, 2>, 2>;

    static constexpr array_t<unsigned int, 2> hermit1 = {0u, 0u};
    static constexpr neighborhood<unsigned int> hermit2 = {hermit1, hermit1};

    grid2() = default;

    DETRAY_HOST
    grid2(vecmem::memory_resource &mr,
          const bare_value m_invalid =
              detray::detail::invalid_value<bare_value>())
        : _axis_p0(mr),
          _axis_p1(mr),
          _data_serialized(&mr),
          _populator(m_invalid) {}

    /** Constructor from axes - copy semantics
     *
     * @param axis_p0 is the axis in the first coordinate
     * @param axis_p1 is the axis in the second coordinate
     *
     **/
    DETRAY_HOST
    grid2(const axis_p0_type &axis_p0, const axis_p1_type &axis_p1,
          vecmem::memory_resource &mr,
          const bare_value m_invalid =
              detray::detail::invalid_value<bare_value>())
        : _axis_p0(axis_p0, mr),
          _axis_p1(axis_p1, mr),
          _data_serialized(&mr),
          _populator(m_invalid) {
        _data_serialized = serialized_storage(_axis_p0.bins() * _axis_p1.bins(),
                                              _populator.init());
    }

    /** Constructor from axes - move semantics
     *
     * @param axis_p0 is the axis in the first coordinate
     * @param axis_p1 is the axis in the second coordinate
     *
     **/
    DETRAY_HOST
    grid2(axis_p0_type &&axis_p0, axis_p1_type &&axis_p1,
          vecmem::memory_resource &mr,
          const bare_value m_invalid =
              detray::detail::invalid_value<bare_value>())
        : _axis_p0(std::move(axis_p0), mr),
          _axis_p1(std::move(axis_p1), mr),
          _data_serialized(&mr),
          _populator(m_invalid) {
        _data_serialized = serialized_storage(_axis_p0.bins() * _axis_p1.bins(),
                                              _populator.init());
    }

    /** Constructor from grid data
     **/
    template <typename grid_view_t>
        requires(!std::is_same_v<grid2, grid_view_t>) &&
                    (!std::is_base_of_v<vecmem::memory_resource, grid_view_t>)
    DETRAY_HOST_DEVICE grid2(const grid_view_t &grid_data,
                             const bare_value m_invalid =
                                 detray::detail::invalid_value<bare_value>())
        : _axis_p0(grid_data._axis_p0_view),
          _axis_p1(grid_data._axis_p1_view),
          _data_serialized(grid_data._data_view),
          _populator(m_invalid) {}

    /** Allow for grid shift, when using a centralized store and indices
     *
     * @param offset is the applied offset shift
     *
     **/
    void shift(const typename populator_type::bare_value &offset) {
        std::ranges::for_each(_data_serialized,
                              [&](auto &ds) { _populator.shift(ds, offset); });
    }

    /** Fill/populate operation
     *
     * @tparam point2_t the 2D local point type
     *
     * @param p2 the point in p2 local frame
     * @param fvalue is a single fill value to be filled
     *
     **/
    template <detray::concepts::point2D point2_t>
    DETRAY_HOST_DEVICE void populate(
        const point2_t &p2, typename populator_type::bare_value &&fvalue) {
        auto sbin = _serializer.template serialize<axis_p0_type, axis_p1_type>(
            _axis_p0, _axis_p1, _axis_p0.bin(p2[0]), _axis_p1.bin(p2[1]));
        _populator(_data_serialized[sbin], std::move(fvalue));
    }

    /** Fill/populate operation - with bin entry
     *
     * @param bin The two-dimensional bin2
     * @param fvalue is a single fill value to be filled
     *
     **/
    DETRAY_HOST_DEVICE
    void populate(unsigned int bin0, unsigned int bin1,
                  typename populator_type::bare_value &&fvalue) {
        auto sbin = _serializer.template serialize<axis_p0_type, axis_p1_type>(
            _axis_p0, _axis_p1, bin0, bin1);
        _populator(_data_serialized[sbin], std::move(fvalue));
    }

    DETRAY_HOST_DEVICE
    void populate(unsigned int gbin,
                  typename populator_type::bare_value &&fvalue) {
        unsigned int bin0 = gbin % _axis_p0.bins();
        unsigned int bin1 = gbin / _axis_p0.bins();
        populate(bin0, bin1, std::move(fvalue));
    }

    /** Return the value of a single bin - with direct bin acess
     *
     * @param bin0 the index of bin 0
     * @param bin1 the index of bin 1
     *
     * @return the const reference to the value in this bin
     **/
    DETRAY_HOST_DEVICE
    typename serialized_storage::const_reference bin(unsigned int bin0,
                                                     unsigned int bin1) const {
        return _data_serialized[_serializer.template serialize<
            axis_p0_type, axis_p1_type>(_axis_p0, _axis_p1, bin0, bin1)];
    }

    DETRAY_HOST_DEVICE
    typename serialized_storage::reference bin(unsigned int bin0,
                                               unsigned int bin1) {
        return _data_serialized.at(
            _serializer.template serialize<axis_p0_type, axis_p1_type>(
                _axis_p0, _axis_p1, bin0, bin1));
    }

    DETRAY_HOST_DEVICE
    typename serialized_storage::const_reference bin(unsigned int gbin) const {
        unsigned int bin0 = gbin % _axis_p0.bins();
        unsigned int bin1 = gbin / _axis_p0.bins();
        return bin(bin0, bin1);
    }

    DETRAY_HOST_DEVICE
    typename serialized_storage::reference bin(unsigned int gbin) {
        unsigned int bin0 = gbin % _axis_p0.bins();
        unsigned int bin1 = gbin / _axis_p0.bins();
        return bin(bin0, bin1);
    }

    /** Return the value of a single bin
     *
     * @param p2 is point in the local frame
     *
     * @return the const reference to the value in this bin
     **/
    template <detray::concepts::point2D point2_t>
        requires(!std::is_scalar_v<point2_t>)
    DETRAY_HOST_DEVICE typename serialized_storage::const_reference bin(
        const point2_t &p2) const {
        return _data_serialized[_serializer.template serialize<axis_p0_type,
                                                               axis_p1_type>(
            _axis_p0, _axis_p1, _axis_p0.bin(p2[0]), _axis_p1.bin(p2[1]))];
    }

    /** Return the value of a single bin - non-const access
     *
     * @param p2 is point in the local frame
     *
     * @return the const reference to the value in this bin
     **/
    template <detray::concepts::point2D point2_t>
        requires(!std::is_scalar_v<point2_t>)
    DETRAY_HOST_DEVICE typename serialized_storage::reference bin(
        const point2_t &p2) {
        return _data_serialized[_serializer.template serialize<axis_p0_type,
                                                               axis_p1_type>(
            _axis_p0, _axis_p1, _axis_p0.bin(p2[0]), _axis_p1.bin(p2[1]))];
    }

    /** Return a zone around a single bin, either with binned or scalar
     *neighborhood
     *
     * The zone is done with a neighborhood around the bin which is defined by
     *p2
     *
     * @param p2 is point in the local frame
     * @param nhood is the binned/scalar neighborhood
     * @param sort is a directive whether to sort or not
     *
     * @return the sequence of values
     **/
    template <typename neighbor_t, detray::concepts::point2D point2_t>
    DETRAY_HOST_DEVICE vector_t<typename populator_type::bare_value> zone_t(
        const point2_t &p2, const neighborhood<neighbor_t> &nhood,
        bool sort) const {
        auto zone0 = _axis_p0.zone(p2[0], nhood[0]);
        auto zone1 = _axis_p1.zone(p2[1], nhood[1]);

        vector_t<typename populator_type::bare_value> zone;

        // Specialization for bare value equal to store value
        if constexpr (std::is_same_v<typename populator_type::bare_value,
                                     typename populator_type::store_value>) {
            unsigned int iz = 0u;
            zone = vector_t<typename populator_type::bare_value>(
                zone0.size() * zone1.size(), {});
            for (const auto z1 : zone1) {
                for (const auto z0 : zone0) {
                    auto sbin =
                        _serializer
                            .template serialize<axis_p0_type, axis_p1_type>(
                                _axis_p0, _axis_p1, z0, z1);
                    zone[iz++] = _data_serialized[sbin];
                }
            }
        } else {
            zone.reserve(10u);
            for (const auto z1 : zone1) {
                for (const auto z0 : zone0) {
                    auto sbin =
                        _serializer
                            .template serialize<axis_p0_type, axis_p1_type>(
                                _axis_p0, _axis_p1, z0, z1);
                    auto bin_data = _data_serialized[sbin];
                    auto bin_content = _populator.sequence(bin_data);
                    zone.insert(zone.end(), bin_content.begin(),
                                bin_content.end());
                }
            }
        }

        if (sort) {
            std::ranges::sort(zone);
        }
        return zone;
    }

    /** Return a zone around a single bin, either with binned neighborhood
     *
     * The zone is done with a neighborhood around the bin which is defined by
     *p2
     *
     * @param p2 is point in the local frame
     * @param nhood is the binned neighborhood
     * @param sort is a directive whether to sort or not
     *
     * @return the sequence of values
     **/
    template <detray::concepts::point2D point2_t>
    DETRAY_HOST_DEVICE vector_t<typename populator_type::bare_value> zone(
        const point2_t &p2, const neighborhood<unsigned int> &nhood = hermit2,
        bool sort = false) const {
        return zone_t<unsigned int>(p2, nhood, sort);
    }

    /** Return a zone around a single bin, either with scalar neighborhood
     *
     * The zone is done with a neighborhood around the bin which is defined by
     *p2
     *
     * @param p2 is point in the local frame
     * @param nhood is the binned neighborhood
     * @param sort is a directive whether to sort or not
     *
     * @return the sequence of values
     **/
    template <detray::concepts::point2D point2_t>
    DETRAY_HOST_DEVICE vector_t<typename populator_type::bare_value> zone(
        const point2_t &p2, const neighborhood<scalar> &nhood,
        bool sort = false) const {
        return zone_t<scalar>(p2, nhood, sort);
    }

    /** Const access to axis p0  */
    DETRAY_HOST_DEVICE
    const axis_p0_type &axis_p0() const { return _axis_p0; }

    /** Non-const access to axis p0  */
    DETRAY_HOST_DEVICE
    axis_p0_type &axis_p0() { return _axis_p0; }

    /** Const access to axis p1 */
    DETRAY_HOST_DEVICE
    const axis_p1_type &axis_p1() const { return _axis_p1; }

    /** Non-const access to axis p1 */
    DETRAY_HOST_DEVICE
    axis_p1_type &axis_p1() { return _axis_p1; }

    /* Copy of axes in a tuple */
    DETRAY_HOST_DEVICE
    tuple_t<axis_p0_type, axis_p1_type> axes() const {
        return std::tie(_axis_p0, _axis_p1);
    }

    /* Get the total number of bins */
    DETRAY_HOST_DEVICE
    unsigned int nbins() const { return _axis_p0.bins() * _axis_p1.bins(); }

    /** Const acess to the serializer */
    DETRAY_HOST_DEVICE
    const serializer_type &serializer() const { return _serializer; }

    /** Const acess to the polulator */
    DETRAY_HOST_DEVICE
    const populator_type &populator() const { return _populator; }

    DETRAY_HOST_DEVICE
    serialized_storage &data() { return _data_serialized; }
    DETRAY_HOST_DEVICE
    const serialized_storage &data() const { return _data_serialized; }

    private:
    axis_p0_type _axis_p0;
    axis_p1_type _axis_p1;
    serialized_storage _data_serialized;
    populator_type _populator;
    serializer_type _serializer;
};

/// Helper function creating a @c vecmem::data::vector_view object
template <typename TYPE, typename ALLOC>
DETRAY_HOST vecmem::data::vector_view<TYPE> get_data(
    std::vector<TYPE, ALLOC> &vec, vecmem::memory_resource &) {

    return vecmem::get_data(vec);
}

/// Helper function creating a @c vecmem::data::vector_view object
template <typename TYPE, typename ALLOC>
DETRAY_HOST vecmem::data::vector_view<const TYPE> get_data(
    const std::vector<TYPE, ALLOC> &vec, vecmem::memory_resource &) {

    return vecmem::get_data(vec);
}

/// Helper function creating a @c vecmem::data::jagged_vector_data object
template <typename TYPE, typename ALLOC1, typename ALLOC2>
DETRAY_HOST vecmem::data::jagged_vector_data<TYPE> get_data(
    std::vector<std::vector<TYPE, ALLOC1>, ALLOC2> &vec,
    vecmem::memory_resource &resource) {

    return vecmem::get_data(vec, &resource);
}

/// Helper function creating a @c vecmem::data::vector_view object
template <typename TYPE, typename ALLOC1, typename ALLOC2>
DETRAY_HOST vecmem::data::jagged_vector_data<const TYPE> get_data(
    const std::vector<std::vector<TYPE, ALLOC1>, ALLOC2> &vec,
    vecmem::memory_resource &resource) {

    return vecmem::get_data(vec, &resource);
}

/** A two-dimensional (non-const) grid view for gpu device usage
 **/
template <typename grid2_t>
struct grid2_view {
    axis_data<typename grid2_t::axis_p0_type, scalar> _axis_p0_view;
    axis_data<typename grid2_t::axis_p1_type, scalar> _axis_p1_view;
    typename grid2_t::populator_type::vector_view_type _data_view;
};

/** A two-dimensional (const) grid view for gpu device usage
 **/
template <typename grid2_t>
struct const_grid2_view {

    /// Declare that a default constructor can/should be generated
    const_grid2_view() = default;
    /// Constructor with the 3 member variables
    DETRAY_HOST_DEVICE
    const_grid2_view(
        const axis_data<typename grid2_t::axis_p0_type, const scalar>
            &axis_p0_view,
        const axis_data<typename grid2_t::axis_p1_type, const scalar>
            &axis_p1_view,
        const typename grid2_t::populator_type::const_vector_view_type
            &data_view)
        : _axis_p0_view(axis_p0_view),
          _axis_p1_view(axis_p1_view),
          _data_view(data_view) {}
    /// Construct a const data object from a non-const one
    DETRAY_HOST_DEVICE
    const_grid2_view(const grid2_view<grid2_t> &parent)
        : _axis_p0_view(parent._axis_p0_view),
          _axis_p1_view(parent._axis_p1_view),
          _data_view(parent._data_view) {}

    axis_data<typename grid2_t::axis_p0_type, const scalar> _axis_p0_view;
    axis_data<typename grid2_t::axis_p1_type, const scalar> _axis_p1_view;
    typename grid2_t::populator_type::const_vector_view_type _data_view;
};

/** A two-dimensional (non-const) grid data for gpu device usage
 **/
template <typename grid2_t>
struct grid2_data : public grid2_view<grid2_t> {

    /** Constructor from grid
     *
     * @param grid is the input grid from host
     * @param resource is the vecmem memory resource
     *
     **/
    grid2_data(grid2_t &grid, vecmem::memory_resource &resource)
        : _axis_p0(grid.axis_p0(), resource),
          _axis_p1(grid.axis_p1(), resource),
          _data(traccc::get_data(grid.data(), resource)) {
        grid2_view<grid2_t>::_axis_p0_view = traccc::get_data(_axis_p0);
        grid2_view<grid2_t>::_axis_p1_view = traccc::get_data(_axis_p1);
        grid2_view<grid2_t>::_data_view = _data;
    }

    typename grid2_t::axis_p0_type _axis_p0;
    typename grid2_t::axis_p1_type _axis_p1;
    typename grid2_t::populator_type::vector_data_type _data;
};

/** A two-dimensional (const) grid data for gpu device usage
 **/
template <typename grid2_t>
struct const_grid2_data : public const_grid2_view<grid2_t> {

    /** Constructor from grid
     *
     * @param grid is the input grid from host
     * @param resource is the vecmem memory resource
     *
     **/
    const_grid2_data(const grid2_t &grid, vecmem::memory_resource &resource)
        : _axis_p0(grid.axis_p0(), resource),
          _axis_p1(grid.axis_p1(), resource),
          _data(traccc::get_data(grid.data(), resource)) {
        const typename grid2_t::axis_p0_type &const_axis_p0 = _axis_p0;
        const_grid2_view<grid2_t>::_axis_p0_view =
            traccc::get_data(const_axis_p0);
        const typename grid2_t::axis_p1_type &const_axis_p1 = _axis_p1;
        const_grid2_view<grid2_t>::_axis_p1_view =
            traccc::get_data(const_axis_p1);
        const_grid2_view<grid2_t>::_data_view = _data;
    }

    typename grid2_t::axis_p0_type _axis_p0;
    typename grid2_t::axis_p1_type _axis_p1;
    typename grid2_t::populator_type::const_vector_data_type _data;
};

/** A two-dimensional grid buffer
 **/
template <typename grid2_t>
struct grid2_buffer : public grid2_view<grid2_t> {

    using populator_type = typename grid2_t::populator_type;
    using axis_p0_type = typename grid2_t::axis_p0_type;
    using axis_p1_type = typename grid2_t::axis_p1_type;

    /** Constructor
     *
     * @param axis_p0 is the first axis
     * @param axis_p1 is the second axis
     * @param capacities is the capacity of vector
     * @param resource is the vecmem memory resource
     * @param host_resurce is the host accessible memory resource
     * @param buffer_type is whether the buffer is resizable or of fixed size
     *
     **/
    grid2_buffer(const axis_p0_type &axis_p0, const axis_p1_type &axis_p1,
                 typename populator_type::buffer_size_type capacities,
                 vecmem::memory_resource &resource,
                 vecmem::memory_resource *host_resource = nullptr,
                 vecmem::data::buffer_type buffer_type =
                     vecmem::data::buffer_type::fixed_size)
        : _axis_p0(axis_p0),
          _axis_p1(axis_p1),
          _buffer(capacities, resource, host_resource, buffer_type) {
        grid2_view<grid2_t>::_axis_p0_view = traccc::get_data(_axis_p0);
        grid2_view<grid2_t>::_axis_p1_view = traccc::get_data(_axis_p1);
        grid2_view<grid2_t>::_data_view = _buffer;
    }

    axis_p0_type _axis_p0;
    axis_p1_type _axis_p1;
    typename populator_type::vector_buffer_type _buffer;
};

/** Get grid2_data from grid and memory resource
 **/
template <template <template <typename...> class, template <typename...> class,
                    template <typename, std::size_t> class, typename, bool,
                    unsigned int>
          class populator_t,
          template <template <typename, std::size_t> class,
                    template <typename...> class>
          class axis_p0_t,
          template <template <typename, std::size_t> class,
                    template <typename...> class>
          class axis_p1_t,
          typename serializer_t, template <typename...> class vector_t,
          template <typename...> class jagged_vector_t,
          template <typename, std::size_t> class array_t,
          template <typename...> class tuple_t, typename value_t = unsigned int,
          bool kSORT = false, unsigned int kDIM = 1>
inline grid2_data<
    grid2<populator_t, axis_p0_t, axis_p1_t, serializer_t, vector_t,
          jagged_vector_t, array_t, tuple_t, value_t, kSORT, kDIM>>
get_data(grid2<populator_t, axis_p0_t, axis_p1_t, serializer_t, vector_t,
               jagged_vector_t, array_t, tuple_t, value_t, kSORT, kDIM> &grid,
         vecmem::memory_resource &resource) {
    return {grid, resource};
}

/** Get const_grid2_data from grid and memory resource
 **/
template <template <template <typename...> class, template <typename...> class,
                    template <typename, std::size_t> class, typename, bool,
                    unsigned int>
          class populator_t,
          template <template <typename, std::size_t> class,
                    template <typename...> class>
          class axis_p0_t,
          template <template <typename, std::size_t> class,
                    template <typename...> class>
          class axis_p1_t,
          typename serializer_t, template <typename...> class vector_t,
          template <typename...> class jagged_vector_t,
          template <typename, std::size_t> class array_t,
          template <typename...> class tuple_t, typename value_t = unsigned int,
          bool kSORT = false, unsigned int kDIM = 1>
inline const_grid2_data<
    grid2<populator_t, axis_p0_t, axis_p1_t, serializer_t, vector_t,
          jagged_vector_t, array_t, tuple_t, value_t, kSORT, kDIM>>
get_data(
    const grid2<populator_t, axis_p0_t, axis_p1_t, serializer_t, vector_t,
                jagged_vector_t, array_t, tuple_t, value_t, kSORT, kDIM> &grid,
    vecmem::memory_resource &resource) {
    return {grid, resource};
}

}  // namespace traccc

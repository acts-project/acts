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

// Detray include(s)
#include <detray/utils/concepts.hpp>

// VecMem include(s).
#include <vecmem/containers/data/vector_view.hpp>
#include <vecmem/containers/device_vector.hpp>
#include <vecmem/memory/memory_resource.hpp>

// System include(s).
#include <algorithm>
#include <cstddef>
#include <type_traits>

namespace traccc {

namespace axis2 {

/** A regular closed axis.
 *
 * The axis is closed, i.e. each underflow bin is mapped to 0
 * and henceforth each overflow bin is mapped to bins-1
 */
template <template <typename, std::size_t> class array_t = std::array,
          template <typename...> class vector_t = vecmem::vector>
struct regular {
    unsigned int n_bins;
    scalar min;
    scalar max;

    static constexpr unsigned int axis_identifier = 0u;

    /** Defualt constructor for dummy axis **/
    DETRAY_HOST_DEVICE
    regular()
        : n_bins(std::numeric_limits<unsigned int>::max()),
          min(static_cast<scalar>(0.)),
          max(static_cast<scalar>(n_bins)) {}

    /** Constructor with vecmem memory resource **/
    DETRAY_HOST
    explicit regular(vecmem::memory_resource & /*resource*/)
        : n_bins(std::numeric_limits<unsigned int>::max()),
          min(static_cast<scalar>(0.)),
          max(static_cast<scalar>(n_bins)) {}

    /** Constructor with vecmem memory resource **/
    DETRAY_HOST
    regular(unsigned int axis_bins, scalar axis_min, scalar axis_max,
            vecmem::memory_resource & /*resource*/)
        : n_bins(axis_bins), min(axis_min), max(axis_max) {}

    /** Constructor with vecmem memory resource **/
    DETRAY_HOST
    regular(const regular &axis, vecmem::memory_resource & /*resource*/)
        : n_bins(axis.n_bins), min(axis.min), max(axis.max) {}

    /** Constructor with axis_data **/
    template <typename axis_data_t>
        requires(!std::is_same_v<regular, axis_data_t>)
    DETRAY_HOST_DEVICE explicit regular(const axis_data_t &axis_data)
        : n_bins(axis_data.n_bins), min(axis_data.min), max(axis_data.max) {}

    /** Return the number of bins */
    DETRAY_HOST_DEVICE
    unsigned int bins() const { return n_bins; }

    /** Access function to a single bin from a value v
     *
     * @param v is the value for the bin search
     *
     * As the axis is closed it @returns a unsigned int type
     **/
    DETRAY_HOST_DEVICE
    unsigned int bin(scalar v) const {
        auto ibin = static_cast<int>((v - min) / (max - min) *
                                     static_cast<scalar>(n_bins));
        if (ibin >= 0 && ibin < static_cast<int>(n_bins)) {
            return static_cast<unsigned int>(ibin);
        } else {
            if (ibin < 0) {
                return 0;
            } else {
                return n_bins - 1u;
            }
        }
    }

    /** Access function to a range with binned neighborhood
     *
     * @param v is the value for the bin search
     * @param nhood is the neighborhood range (in #bins)
     *
     * As the axis is closed it @returns an index range
     **/
    DETRAY_HOST_DEVICE
    std::array<unsigned int, 2> range(scalar v,
                                      const array_t<unsigned int, 2> &nhood = {
                                          0u, 0u}) const {

        auto ibin = static_cast<int>((v - min) / (max - min) *
                                     static_cast<scalar>(n_bins));
        int ibinmin = ibin - static_cast<int>(nhood[0]);
        int ibinmax = ibin + static_cast<int>(nhood[1]);
        unsigned int min_bin =
            (ibinmin >= 0) ? static_cast<unsigned int>(ibinmin) : 0u;
        unsigned int max_bin = (ibinmax < static_cast<int>(n_bins))
                                   ? static_cast<unsigned int>(ibinmax)
                                   : n_bins - 1u;
        return {min_bin, max_bin};
    }

    /** Access function to a range with scalar neighborhood
     *
     * @param v is the value for the bin search
     * @param nhood is the neighborhood range (scalar)
     *
     * As the axis is closed it @returns an index range
     **/
    DETRAY_HOST_DEVICE
    std::array<unsigned int, 2> range(scalar v,
                                      const array_t<scalar, 2> &nhood) const {
        auto nbin = static_cast<int>((v - nhood[0] - min) / (max - min) *
                                     static_cast<scalar>(n_bins));
        auto pbin = static_cast<int>((v + nhood[1] - min) / (max - min) *
                                     static_cast<scalar>(n_bins));
        unsigned int min_bin =
            (nbin >= 0) ? static_cast<unsigned int>(nbin) : 0u;
        unsigned int max_bin = (pbin < static_cast<int>(n_bins))
                                   ? static_cast<unsigned int>(pbin)
                                   : n_bins - 1;
        return {min_bin, max_bin};
    }

    /** Access function to a zone with binned neighborhood
     *
     *
     * @tparam neighbor_t is the neighborhood size
     *
     * @param v is the value for the bin search
     * @param nhood is the neighborhood range (in #bins)
     *
     * As the axis is closed it @returns an index sequence
     **/
    template <typename neighbor_t>
    DETRAY_HOST_DEVICE vecmem::vector<unsigned int> zone_t(
        scalar v, const array_t<neighbor_t, 2> &nhood) const {
        std::array<unsigned int, 2> nh_range = range(v, nhood);
        vecmem::vector<unsigned int> sequence(
            static_cast<vecmem::vector<unsigned int>::size_type>(
                nh_range[1] - nh_range[0] + 1u),
            nh_range[0]);
        unsigned int m = 0u;
        std::ranges::for_each(sequence, [&](auto &n) { n += m++; });
        return sequence;
    }

    /** Access function to a zone with binned neighborhood
     *
     * @param v is the value for the bin search
     * @param nhood is the neighborhood range (#bins)
     *
     * As the axis is closed it @returns an index sequence
     **/
    DETRAY_HOST_DEVICE
    vecmem::vector<unsigned int> zone(
        scalar v, const array_t<unsigned int, 2> &nhood) const {
        return zone_t<unsigned int>(v, nhood);
    }

    /** Access function to a zone with scalar neighborhood
     *
     * @param v is the value for the bin search
     * @param nhood is the neighborhood range (scalar)
     *
     * As the axis is closed it @returns an index sequence
     **/
    DETRAY_HOST_DEVICE
    vecmem::vector<unsigned int> zone(scalar v,
                                      const array_t<scalar, 2> &nhood) const {
        return zone_t<scalar>(v, nhood);
    }

    /** @return the bin boundaries for a given @param ibin */
    DETRAY_HOST_DEVICE
    array_t<scalar, 2> borders(unsigned int ibin) const {
        scalar step = (max - min) / static_cast<scalar>(n_bins);
        return {min + static_cast<scalar>(ibin) * step,
                min + static_cast<scalar>(ibin + 1u) * step};
    }

    /** @return the values of the borders */
    DETRAY_HOST_DEVICE
    vector_t<scalar> all_borders() const {
        vector_t<scalar> borders;
        borders.reserve(n_bins + 1u);
        scalar step = (max - min) / static_cast<scalar>(n_bins);
        for (unsigned int ib = 0u; ib < n_bins + 1u; ++ib) {
            borders.push_back(min + static_cast<scalar>(ib) * step);
        }
        return borders;
    }

    /** @return the axis span [min, max) */
    DETRAY_HOST_DEVICE
    array_t<scalar, 2> span() const { return {min, max}; }
};

/** A regular circular axis.
 *
 * The axis is circular, i.e. the underflow bins map into the circular sequence
 */
template <template <typename, std::size_t> class array_t = std::array,
          template <typename...> class vector_t = vecmem::vector>
struct circular {

    unsigned int n_bins;
    scalar min;
    scalar max;

    static constexpr unsigned int axis_identifier = 1u;

    /** Defualt constructor for dummy axis **/
    DETRAY_HOST_DEVICE
    circular()
        : n_bins(std::numeric_limits<unsigned int>::max()),
          min(static_cast<scalar>(0.)),
          max(static_cast<scalar>(n_bins)) {}

    /** Constructor with vecmem memory resource **/
    DETRAY_HOST
    explicit circular(vecmem::memory_resource & /*resource*/)
        : n_bins(std::numeric_limits<unsigned int>::max()),
          min(static_cast<scalar>(0.)),
          max(static_cast<scalar>(n_bins)) {}

    /** Constructor with vecmem memory resource **/
    DETRAY_HOST
    circular(unsigned int axis_bins, scalar axis_min, scalar axis_max,
             vecmem::memory_resource & /*resource*/)
        : n_bins(axis_bins), min(axis_min), max(axis_max) {}

    /** Constructor with vecmem memory resource **/
    DETRAY_HOST
    circular(const circular &axis, vecmem::memory_resource & /*resource*/)
        : n_bins(axis.n_bins), min(axis.min), max(axis.max) {}

    /** Constructor with axis_data **/
    template <typename axis_data_t>
        requires(!std::is_same_v<circular, axis_data_t>)
    DETRAY_HOST_DEVICE explicit circular(const axis_data_t &axis_data)
        : n_bins(axis_data.n_bins), min(axis_data.min), max(axis_data.max) {}

    /** Return the number of bins */
    DETRAY_HOST_DEVICE
    unsigned int bins() const { return n_bins; }

    /** Access function to a single bin from a value v
     *
     * @param v is the value for the bin search
     *
     * As the axis is closed it @returns an index
     **/
    DETRAY_HOST_DEVICE
    unsigned int bin(scalar v) const {
        auto ibin = static_cast<int>((v - min) / (max - min) *
                                     static_cast<scalar>(n_bins));
        if (ibin >= 0 && ibin < static_cast<int>(n_bins)) {
            return static_cast<unsigned int>(ibin);
        } else {
            if (ibin < 0) {
                return n_bins + static_cast<unsigned int>(ibin);
            } else {
                return static_cast<unsigned int>(ibin) - n_bins;
            }
        }
    }

    /** Access function to a range with binned neighborhood
     *
     * @param v is the value for the bin search
     * @param nhood is the  neighborhood range (in #bins)
     *
     * As the axis is circular it @returns an index range
     **/
    DETRAY_HOST_DEVICE
    std::array<unsigned int, 2> range(scalar v,
                                      const array_t<unsigned int, 2> nhood = {
                                          0u, 0u}) const {
        unsigned int gbin = bin(v);
        unsigned int min_bin = remap(gbin, -static_cast<int>(nhood[0]));
        unsigned int max_bin = remap(gbin, static_cast<int>(nhood[1]));
        return {min_bin, max_bin};
    }

    /** Access function to a range with scalar neighborhood
     *
     * @param v is the value for the bin search
     * @param nhood is the neighborhood range (scalar)
     *
     * As the axis is circular it @returns an index range
     **/
    DETRAY_HOST_DEVICE
    std::array<unsigned int, 2> range(scalar v,
                                      const array_t<scalar, 2> &nhood) const {
        unsigned int nbin = bin(v - nhood[0]);
        unsigned int pbin = bin(v + nhood[1]);
        return {nbin, pbin};
    }

    /** Access function to a zone with binned/scalar neighborhood
     *
     * @tparam neighbor_t is the neighborhood size
     *
     * @param v is the value for the bin search
     * @param nhood is the neighborhood range (in #bins/scalar)
     *
     * As the axis is closed it @returns an index sequence
     **/
    template <typename neighbor_t>
    DETRAY_HOST_DEVICE vecmem::vector<unsigned int> zone_t(
        scalar v, const array_t<neighbor_t, 2> &nhood) const {
        std::array<unsigned int, 2> nh_range = range(v, nhood);
        if (nh_range[0] < nh_range[1]) {
            vecmem::vector<unsigned int> sequence(
                static_cast<vecmem::vector<unsigned int>::size_type>(
                    nh_range[1] - nh_range[0] + 1u),
                nh_range[0]);
            unsigned int m = 0;
            std::ranges::for_each(sequence, [&](auto &n) { n += m++; });
            return sequence;
        }
        unsigned int vl = n_bins - nh_range[0] + nh_range[1] + 1u;
        unsigned int mi = 0;
        unsigned int mo = 0;
        vecmem::vector<unsigned int> sequence(
            static_cast<vecmem::vector<unsigned int>::size_type>(vl),
            nh_range[0]);
        std::ranges::for_each(sequence, [&](auto &n) {
            n += mi++;
            if (n > n_bins - 1u) {
                n = mo++;
            }
        });
        return sequence;
    }

    /** Access function to a zone with binned neighborhood
     *
     * @param v is the value for the bin search
     * @param nhood is the neighborhood range (#bins)
     *
     * As the axis is closed it @returns an index sequence
     **/
    DETRAY_HOST_DEVICE
    vecmem::vector<unsigned int> zone(
        scalar v, const array_t<unsigned int, 2> &nhood) const {
        return zone_t<unsigned int>(v, nhood);
    }

    /** Access function to a zone with scalar neighborhood
     *
     * @param v is the value for the bin search
     * @param nhood is the neighborhood range (scalar)
     *
     * As the axis is closed it @returns an index sequence
     **/
    DETRAY_HOST_DEVICE
    vecmem::vector<unsigned int> zone(scalar v,
                                      const array_t<scalar, 2> &nhood) const {
        return zone_t<scalar>(v, nhood);
    }

    /** Helper function to remap onto a circular range
     *
     * @param ibin is the optional binning value
     * @param shood is the sided neighbour hood
     *
     * @return an index, remapped bin
     **/
    DETRAY_HOST_DEVICE
    unsigned int remap(unsigned int ibin, int shood) const {
        int opt_bin = static_cast<int>(ibin) + shood;
        if (opt_bin >= 0 && opt_bin < static_cast<int>(n_bins)) {
            return static_cast<unsigned int>(opt_bin);
        }
        if (opt_bin < 0) {
            return static_cast<unsigned int>(static_cast<int>(n_bins) +
                                             opt_bin);
        }
        return static_cast<unsigned int>(opt_bin) - n_bins;
    }

    /** @return the bin boundaries for a given @param ibin */
    DETRAY_HOST_DEVICE
    array_t<scalar, 2> borders(unsigned int ibin) const {
        scalar step = (max - min) / n_bins;
        return {min + static_cast<scalar>(ibin) * step,
                min + static_cast<scalar>(ibin + 1u) * step};
    }

    /** @return the values of the borders */
    DETRAY_HOST_DEVICE
    vector_t<scalar> all_borders() const {
        vector_t<scalar> borders;
        borders.reserve(n_bins + 1u);
        scalar step = (max - min) / static_cast<scalar>(n_bins);
        for (unsigned int ib = 0u; ib < n_bins + 1u; ++ib) {
            borders.push_back(min + static_cast<scalar>(ib) * step);
        }
        return borders;
    }

    /** @return the range  */
    DETRAY_HOST_DEVICE
    array_t<scalar, 2> span() const { return {min, max}; }
};

/** An iregular circular axis.
 *
 * The axis is closed, i.e. the underflow is mapped into the first,
 * the overflow is mapped into the last.
 */
template <template <typename, std::size_t> class array_t = std::array,
          template <typename...> class vector_t = vecmem::vector>
struct irregular {

    /* dummy bin size, min and max */
    unsigned int n_bins;
    scalar min;
    scalar max;

    vector_t<scalar> boundaries;

    static constexpr unsigned int axis_identifier = 2u;

    /** Defualt constructor for dummy axis **/
    DETRAY_HOST_DEVICE
    irregular()
        : n_bins(std::numeric_limits<unsigned int>::max()),
          min(0.f),
          max(static_cast<scalar>(n_bins)),
          boundaries({}) {}

    /** Constructor with vecmem memory resource **/
    DETRAY_HOST
    explicit irregular(vecmem::memory_resource &resource)
        : n_bins(std::numeric_limits<unsigned int>::max()),
          min(0.f),
          max(static_cast<scalar>(n_bins)),
          boundaries(&resource) {}

    /** Constructor with vecmem memory resource - rvalue **/
    DETRAY_HOST irregular(vector_t<scalar> &&bins,
                          vecmem::memory_resource &resource)
        : n_bins(static_cast<unsigned int>(bins.size()) - 1u),
          min(bins[0]),
          max(bins[n_bins]),
          boundaries(bins, &resource) {}

    /** Constructor with vecmem memory resource - lvalue **/
    DETRAY_HOST irregular(const vector_t<scalar> &bins,
                          vecmem::memory_resource &resource)
        : n_bins(static_cast<unsigned int>(bins.size()) - 1u),
          min(bins[0]),
          max(bins[n_bins]),
          boundaries(bins, &resource) {}

    /** Constructor with vecmem memory resource **/
    DETRAY_HOST
    irregular(const irregular &axis, vecmem::memory_resource &resource)
        : n_bins(axis.n_bins),
          min(axis.min),
          max(axis.max),
          boundaries(axis.boundaries, &resource) {}

    /** Constructor with axis_data **/
    template <typename axis_data_t>
        requires(!std::is_same_v<irregular, axis_data_t>)
    DETRAY_HOST_DEVICE explicit irregular(const axis_data_t &axis_data)
        : n_bins(axis_data.n_bins),
          min(axis_data.min),
          max(axis_data.max),
          boundaries(axis_data.boundaries) {}

    /** Return the number of bins */
    DETRAY_HOST_DEVICE
    unsigned int bins() const {
        return static_cast<unsigned int>(boundaries.size() - 1u);
    }

    /** Access function to a single bin from a value v
     *
     * @param v is the value for the bin search
     *
     * As the axis is closed it @returns a unsigned int type
     **/
    DETRAY_HOST_DEVICE
    unsigned int bin(scalar v) const {
        auto ibin = static_cast<int>(std::ranges::lower_bound(boundaries, v) -
                                     boundaries.begin());
        if (ibin > 0 && ibin < static_cast<int>(boundaries.size())) {
            return static_cast<unsigned int>(--ibin);
        } else {
            if (ibin == 0) {
                return static_cast<unsigned int>(ibin);
            } else {
                return static_cast<unsigned int>(boundaries.size() - 2u);
            }
        }
    }

    /** Access function to a range with binned neighborhood
     *
     * @param v is the value for the bin search
     * @param nhood is the neighborhood range (#bins)
     *
     * As the axis is closed it @returns an index range
     **/
    DETRAY_HOST_DEVICE
    std::array<unsigned int, 2> range(scalar v,
                                      const array_t<unsigned int, 2> &nhood = {
                                          0u, 0u}) const {

        unsigned int ibin = bin(v);
        int bins = static_cast<int>(boundaries.size()) - 1;
        int ibinmin = static_cast<int>(ibin) - static_cast<int>(nhood[0]);
        int ibinmax = static_cast<int>(ibin + nhood[1]);
        auto min_bin = (ibinmin >= 0) ? static_cast<unsigned int>(ibinmin) : 0u;
        auto max_bin = (ibinmax < bins) ? static_cast<unsigned int>(ibinmax)
                                        : static_cast<unsigned int>(bins - 1);
        return {min_bin, max_bin};
    }

    /** Access function to a range with scalar neighborhood
     *
     * @param v is the value for the bin search
     * @param nhood is the neighborhood range (scalar)
     *
     * As the axis is closed it @returns an index range
     **/
    DETRAY_HOST_DEVICE
    std::array<unsigned int, 2> range(scalar v,
                                      const array_t<scalar, 2> &nhood) const {
        unsigned int nbin = bin(v - nhood[0]);
        unsigned int pbin = bin(v + nhood[1]);
        return {nbin, pbin};
    }

    /** Access function to a zone with binned/scalar neighborhood
     *
     * @tparam neighbor_t is the neighborhood type
     *
     * @param v is the value for the bin search
     * @param nhood is the neighborhood range (binned/scalar)
     *
     * As the axis is closed it @returns an index sequence
     **/
    template <typename neighbor_t>
    DETRAY_HOST_DEVICE vecmem::vector<unsigned int> zone_t(
        scalar v, const array_t<neighbor_t, 2> nhood) const {
        std::array<unsigned int, 2> nh_range = range(v, nhood);
        vecmem::vector<unsigned int> sequence(
            static_cast<vecmem::vector<unsigned int>::size_type>(
                nh_range[1] - nh_range[0] + 1u),
            nh_range[0]);
        unsigned int m = 0u;
        std::ranges::for_each(sequence, [&](auto &n) { n += m++; });
        return sequence;
    }

    /** Access function to a zone with binned neighborhood
     *
     * @param v is the value for the bin search
     * @param nhood is the neighborhood range (#bins)
     *
     * As the axis is closed it @returns an index sequence
     **/
    DETRAY_HOST_DEVICE
    vecmem::vector<unsigned int> zone(scalar v,
                                      const array_t<unsigned int, 2> &nhood = {
                                          0u, 0u}) const {
        return zone_t<unsigned int>(v, nhood);
    }

    /** Access function to a zone with scalar neighborhood
     *
     * @param v is the value for the bin search
     * @param nhood is the neighborhood range (scalar)
     *
     * As the axis is closed it @returns an index sequence
     **/
    DETRAY_HOST_DEVICE
    vecmem::vector<unsigned int> zone(scalar v,
                                      const array_t<scalar, 2> &nhood) const {
        return zone_t<scalar>(v, nhood);
    }

    DETRAY_HOST_DEVICE
    /** @return the bin boundaries for a given @param ibin */
    array_t<scalar, 2> borders(unsigned int ibin) const {
        return {boundaries[ibin], boundaries[ibin + 1u]};
    }

    /** @return the values of the borders of all bins */
    DETRAY_HOST
    vector_t<scalar> all_borders() const { return boundaries; }

    /** @return the range  */
    DETRAY_HOST_DEVICE
    array_t<scalar, 2> span() const {
        return {boundaries[0], boundaries[boundaries.size() - 1u]};
    }
};

}  // namespace axis2

/**
 * static implementation of axis data for device
 */
template <typename axis_t, detray::concepts::scalar scalar_t,
          typename Enable = void>
struct axis_data;

template <typename axis_t, detray::concepts::scalar scalar_t>
    requires(axis_t::axis_identifier == 0u) || (axis_t::axis_identifier == 1u)
struct axis_data<axis_t, scalar_t> {

    /// Declare that a default constructor can/should be generated
    axis_data() = default;
    /// Constructor with the 3 member values
    DETRAY_HOST_DEVICE
    axis_data(unsigned int _n_bins, std::remove_cv_t<scalar_t> _min,
              std::remove_cv_t<scalar_t> _max)
        : n_bins(_n_bins), min(_min), max(_max) {}
    /// Construct a const data object from a non-const one
    template <typename other_scalar_t>
        requires detray::concepts::same_as_no_const<scalar_t, other_scalar_t>
    DETRAY_HOST_DEVICE explicit axis_data(
        const axis_data<axis_t, other_scalar_t, void> &parent)
        : n_bins(parent.n_bins), min(parent.min), max(parent.max) {}

    unsigned int n_bins;
    std::remove_cv_t<scalar_t> min;
    std::remove_cv_t<scalar_t> max;
};

template <typename axis_t, detray::concepts::scalar scalar_t>
    requires(axis_t::axis_identifier == 2)
struct axis_data<axis_t, scalar_t> {

    /// Declare that a default constructor can/should be generated
    axis_data() = default;
    /// Constructor with the 4 member values
    DETRAY_HOST_DEVICE
    axis_data(unsigned int _n_bins, std::remove_cv_t<scalar_t> _min,
              std::remove_cv_t<scalar_t> _max,
              const vecmem::data::vector_view<scalar_t> &_boundaries)
        : n_bins(_n_bins), min(_min), max(_max), boundaries(_boundaries) {}
    /// Construct a const data object from a non-const one
    template <typename other_scalar_t>
        requires detray::concepts::same_as_no_const<scalar_t, other_scalar_t>
    DETRAY_HOST_DEVICE explicit axis_data(
        const axis_data<axis_t, other_scalar_t, void> &parent)
        : n_bins(parent.n_bins),
          min(parent.min),
          max(parent.max),
          boundaries(parent.boundaries) {}

    unsigned int n_bins;
    std::remove_cv_t<scalar_t> min;
    std::remove_cv_t<scalar_t> max;
    vecmem::data::vector_view<scalar_t> boundaries;
};

/**
 * standalone function to get axis_data (non-const)
 */
template <template <template <typename, std::size_t> class,
                    template <typename...> class>
          class axis_t,
          template <typename, std::size_t> class array_t,
          template <typename...> class vector_t>
    requires(axis_t<array_t, vector_t>::axis_identifier == 0u) ||
            (axis_t<array_t, vector_t>::axis_identifier == 1u)
inline axis_data<axis_t<array_t, vector_t>, scalar> get_data(
    axis_t<array_t, vector_t> &axis) {

    axis_data<axis_t<array_t, vector_t>, scalar> result{axis.n_bins, axis.min,
                                                        axis.max};
    return result;
}

/**
 * standalone function to get axis_data (non-const)
 */
template <template <template <typename, std::size_t> class,
                    template <typename...> class>
          class axis_t,
          template <typename, std::size_t> class array_t,
          template <typename...> class vector_t>
    requires(axis_t<array_t, vector_t>::axis_identifier == 2u)
inline axis_data<axis_t<array_t, vector_t>, scalar> get_data(
    axis_t<array_t, vector_t> &axis) {

    axis_data<axis_t<array_t, vector_t>, scalar> result{
        axis.n_bins, axis.min, axis.max, vecmem::get_data(axis.boundaries)};
    return result;
}

/**
 * standalone function to get axis_data (const)
 */
template <template <template <typename, std::size_t> class,
                    template <typename...> class>
          class axis_t,
          template <typename, std::size_t> class array_t,
          template <typename...> class vector_t>
    requires(axis_t<array_t, vector_t>::axis_identifier == 0u) ||
            (axis_t<array_t, vector_t>::axis_identifier == 1u)
inline axis_data<axis_t<array_t, vector_t>, const scalar> get_data(
    const axis_t<array_t, vector_t> &axis) {

    axis_data<axis_t<array_t, vector_t>, const scalar> result{
        axis.n_bins, axis.min, axis.max};
    return result;
}

/**
 * standalone function to get axis_data (const)
 */
template <template <template <typename, std::size_t> class,
                    template <typename...> class>
          class axis_t,
          template <typename, std::size_t> class array_t,
          template <typename...> class vector_t>
    requires(axis_t<array_t, vector_t>::axis_identifier == 2u)
inline axis_data<axis_t<array_t, vector_t>, const scalar> get_data(
    const axis_t<array_t, vector_t> &axis) {

    axis_data<axis_t<array_t, vector_t>, const scalar> result{
        axis.n_bins, axis.min, axis.max, vecmem::get_data(axis.boundaries)};
    return result;
}

}  // namespace traccc

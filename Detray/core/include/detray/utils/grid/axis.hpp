// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// Project include(s).
#include "detray/core/detail/container_buffers.hpp"
#include "detray/core/detail/container_views.hpp"
#include "detray/definitions/detail/qualifiers.hpp"
#include "detray/definitions/grid_axis.hpp"
#include "detray/definitions/indexing.hpp"
#include "detray/utils/grid/concepts.hpp"
#include "detray/utils/grid/detail/axis_binning.hpp"
#include "detray/utils/grid/detail/axis_bounds.hpp"
#include "detray/utils/ranges.hpp"
#include "detray/utils/type_list.hpp"
#include "detray/utils/type_registry.hpp"
#include "detray/utils/type_traits.hpp"

// VecMem include(s).
#include <vecmem/memory/memory_resource.hpp>

// System include(s).
#include <cstddef>
#include <type_traits>

namespace detray::axis {

/// @brief Multi-bin: contains bin indices from multiple axes
template <std::size_t DIM, typename index_t = dindex>
struct multi_bin : public dmulti_index<index_t, DIM> {
  using base_t = dmulti_index<index_t, DIM>;
  using base_t::base_t;
};

/// @brief Helper to tie two bin indices to a range.
/// @note Cannot use dindex_range for signed integer bin indices.
using bin_range = darray<int, 2>;

/// @brief Multi-bin-range: contains bin index ranges from multiple axes
template <std::size_t DIM>
struct multi_bin_range : public dmulti_index<bin_range, DIM> {
  using base_t = dmulti_index<bin_range, DIM>;
  using base_t::base_t;
};

/// @brief A single axis.
///
/// An axis ties bounds and binning behaviour with the bin edges storage.
/// The bounds determine how bin indices are mapped at the over-
/// and underflow bins. The type of binning determines whether the axis has
/// regular or irregular binning. The bin edges needed to find bin indices
/// are not owned by the axis, but are passed to the binning type.
template <typename bounds_t, typename binning_t>
struct single_axis {
  /// Make axis bounds accessible
  using bounds_type = bounds_t;
  /// Make axis binning type accessible
  using binning_type = binning_t;

  using scalar_type = typename binning_type::scalar_type;

  /// Extract container types
  /// @{
  using container_types = typename binning_type::container_types;
  template <typename T>
  using vector_type = typename binning_type::template vector_type<T>;
  /// @}

  /// Defines the geometrical bounds of the axis as a service:
  /// open, closed or circular
  DETRAY_NO_UNIQUE_ADDRESS bounds_type m_bounds{};
  /// Defines the binning on the axis as a service: regular vs irregular
  binning_type m_binning{};

  /// Parameterized constructor - empty binning
  constexpr single_axis() = default;

  /// Parametrized constructor that builds the binning scheme
  template <typename... Args>
  DETRAY_HOST_DEVICE single_axis(const dsized_index_range &indx_range,
                                 const vector_type<scalar_type> *edges)
      : m_binning(indx_range, edges) {}

  /// Equality operator
  ///
  /// @param rhs is the right-hand side of the comparison
  ///
  /// @returns whether the two axes are equal
  constexpr bool operator==(const single_axis &rhs) const = default;

  /// @returns the axis label, i.e. x, y, z, r or phi axis.
  DETRAY_HOST_DEVICE
  constexpr auto label() const -> axis::label { return bounds_type::label; }

  /// @returns the type of bounds of the axis, i.e. closed, open or circular.
  DETRAY_HOST_DEVICE
  constexpr auto bounds() const -> axis::bounds { return bounds_type::type; }

  /// @returns the type of binning of the axis, i.e. regular or irregular.
  DETRAY_HOST_DEVICE
  constexpr auto binning() const -> axis::binning { return binning_type::type; }

  /// @returns the total number of bins
  DETRAY_HOST_DEVICE
  constexpr dindex nbins() const {
    // The open axis boundary has extra over- and underflow bins that are
    // automatically added beyond the axis span
    if constexpr (bounds_type::type == axis::bounds::e_open) {
      return m_binning.nbins() + 2u;
    } else {
      return m_binning.nbins();
    }
  }

  /// @returns the width of a bin
  template <typename... Args>
  DETRAY_HOST_DEVICE constexpr scalar_type bin_width(Args &&...args) const {
    return m_binning.bin_width(std::forward<Args>(args)...);
  }

  /// Given a value on the axis, find the correct bin.
  ///
  /// @note This includes bin index wrapping for circular axis.
  ///
  /// @param v is the value for the bin search
  ///
  /// @returns the bin index.
  DETRAY_HOST_DEVICE
  dindex bin(const scalar_type v) const {
    int b{m_bounds.map(m_binning.bin(v), m_binning.nbins())};

    if constexpr (bounds_type::type == axis::bounds::e_circular) {
      b = m_bounds.wrap(b, m_binning.nbins());
    }

    return static_cast<dindex>(b);
  }

  /// Given a value on the axis and a neighborhood, find the correct bin range
  ///
  /// @note (!) The circular axis index wrap-around happens in a separate
  ///           step so that index sequences are well defined
  ///
  /// @param v is the value for the bin search
  /// @param nhood is the neighborhood range (in #bins or value interval)
  ///
  /// @returns a dindex_range around the bin index.
  template <typename neighbor_t>
  DETRAY_HOST_DEVICE bin_range range(const scalar_type v,
                                     const darray<neighbor_t, 2> &nhood) const {
    return m_bounds.map(m_binning.range(v, nhood), m_binning.nbins());
  }

  /// @returns the bin edges for a given @param ibin .
  DETRAY_HOST_DEVICE
  darray<scalar_type, 2> bin_edges(const dindex ibin) const {
    return m_binning.bin_edges(ibin);
  }

  /// @returns the values of the bin edges. Is a succession of lower edges.
  DETRAY_HOST_DEVICE
  vector_type<scalar_type> bin_edges() const { return m_binning.bin_edges(); }

  /// @returns the axis span [min, max).
  DETRAY_HOST_DEVICE
  darray<scalar_type, 2> span() const { return m_binning.span(); }

  /// @returns the axis span [min, max).
  DETRAY_HOST_DEVICE
  scalar_type min() const { return m_binning.span()[0]; }

  /// @returns the axis span [min, max).
  DETRAY_HOST_DEVICE
  scalar_type max() const { return m_binning.span()[1]; }

  /// @returns a string stream that prints the single axis details
  DETRAY_HOST
  friend std::ostream &operator<<(std::ostream &os, const single_axis &ax) {
    os << "label: " << ax.label() << std::endl;
    os << "bounds: " << ax.bounds() << std::endl;
    os << "binning: " << ax.binning() << std::endl;
    os << "n-bins: " << ax.nbins() << std::endl;
    os << "span: [" << ax.min() << ", " << ax.max() << "]";

    return os;
  }
};

/// @brief An N-dimensional collection of single axes.
///
/// Given a point in the grid local coordinate system, which is spanned by the
/// axes in this multi-axes type, the corresponding bin multi-index or
/// multi-index range is returned.
///
/// @note can be owning the data (as member of a standalone grid) or can be
/// non-owning if the grid is part of a larger collection.
template <bool ownership, typename local_frame_t, concepts::axis... axis_ts>
class multi_axis {
  /// Match an axis to its label at compile time
  using axis_reg = types::registry<axis::label, axis_ts...>;

 public:
  /// Dimension of the local coordinate system that is spanned by the axes
  static constexpr dindex dim = sizeof...(axis_ts);
  static constexpr bool is_owning = ownership;

  /// binnings and axis bounds
  /// @{
  using binnings = types::list<typename axis_ts::binning_type...>;
  using bounds = types::list<typename axis_ts::bounds_type...>;
  using loc_bin_index = axis::multi_bin<dim>;
  /// @}

  /// Projection onto local coordinate system that is spanned by the axes
  using local_frame_type = local_frame_t;
  using algebra_type = typename local_frame_type::algebra_type;
  using point_type = typename local_frame_type::loc_point;

  using scalar_type = typename detray::detail::first_t<axis_ts...>::scalar_type;

  /// Extract container types
  /// @{
  using container_types =
      typename detray::detail::first_t<axis_ts...>::container_types;
  template <typename T>
  using vector_type = typename container_types::template vector_type<T>;
  /// @}

 private:
  /// Owning and non-owning range of edge offsets
  using edge_offset_range_t = std::conditional_t<
      is_owning, vector_type<dsized_index_range>,
      detray::ranges::subrange<const vector_type<dsized_index_range>>>;
  /// Owning and non-owning range of bin edges
  using edge_range_t = std::conditional_t<is_owning, vector_type<scalar_type>,
                                          const vector_type<scalar_type> *>;

 public:
  /// Axes boundary/bin edges storage
  /// @{
  using edge_offset_container_type = vector_type<dsized_index_range>;
  using edges_container_type = vector_type<scalar_type>;
  /// @}

  /// Vecmem based view type
  using view_type =
      dmulti_view<dvector_view<dsized_index_range>, dvector_view<scalar_type>>;
  using const_view_type = dmulti_view<dvector_view<const dsized_index_range>,
                                      dvector_view<const scalar_type>>;
  /// Vecmem based buffer type
  using buffer_type = dmulti_buffer<dvector_buffer<dsized_index_range>,
                                    dvector_buffer<scalar_type>>;

  /// Find the corresponding (non-)owning type
  template <bool owning>
  using type = multi_axis<owning, local_frame_t, axis_ts...>;

  /// Default constructor
  constexpr multi_axis() = default;

  /// Construct containers using a specific memory resources
  template <bool owner = is_owning>
    requires owner
  DETRAY_HOST explicit multi_axis(vecmem::memory_resource &resource)
      : m_edge_offsets(&resource), m_edges(&resource) {}

  /// Construct from containers - move
  template <bool owner = is_owning>
    requires owner
  DETRAY_HOST_DEVICE multi_axis(edge_offset_range_t &&edge_offsets,
                                edge_range_t &&edges)
      : m_edge_offsets(std::move(edge_offsets)), m_edges(std::move(edges)) {}

  /// Construct from containers that are not owned by this class
  ///
  /// @param edge_offsets offsets into the global edge container
  /// @param edges the global edge container
  /// @param offset offset into the global edge offset container
  template <bool owner = is_owning>
    requires(!owner)
  DETRAY_HOST_DEVICE multi_axis(
      const vector_type<dsized_index_range> &edge_offsets,
      const vector_type<scalar_type> &edges, const unsigned int offset = 0)
      : m_edge_offsets(edge_offsets, dsized_index_range{offset, offset + dim}),
        m_edges(&edges) {}

  /// Construct containers from vecmem based view type
  ///
  /// @param view vecmem view on the axes data
  template <concepts::device_view view_t>
  DETRAY_HOST_DEVICE explicit multi_axis(const view_t &view)
      : m_edge_offsets(detray::detail::get<0>(view.m_view)),
        m_edges(detray::detail::get<1>(view.m_view)) {}

  /// @returns access to the underlying bin edge offset storage - const
  DETRAY_HOST_DEVICE
  constexpr auto bin_edge_offsets() const -> const edge_offset_range_t & {
    return m_edge_offsets;
  }

  /// @returns access to the underlying bin edge storage - const
  DETRAY_HOST_DEVICE
  constexpr auto bin_edges() const -> const vector_type<scalar_type> & {
    if constexpr (is_owning) {
      return m_edges;
    } else {
      return *m_edges;
    }
  }

  /// Build an axis object in place.
  /// @{
  /// @tparam I the position of the axis in the parameter pack. Also
  ///               determines which axes data are used to build the instance.
  /// @returns an axis object, corresponding to the index.
  template <std::size_t I>
  DETRAY_HOST_DEVICE types::get<axis_reg, types::id_cast<axis_reg, I>>
  get_axis() const {
    if constexpr (std::same_as<edge_offset_range_t,
                               vecmem::vector<dsized_index_range>>) {
#if defined(__CUDACC__)
      // Otherwise, a warning is triggered with gcc 11.4 and nvcc 12.4
      DETRAY_VERBOSE_DEVICE(
          "The host container types must not be called in device code");
      assert(false);
      return {dsized_index_range{}, &bin_edges()};
#else
      return {m_edge_offsets[I], &bin_edges()};
#endif
    } else {
      return {m_edge_offsets[I], &bin_edges()};
    }
  }

  /// @tparam L label of the axis.
  /// @returns an axis object, corresponding to the label.
  template <axis::label L>
  DETRAY_HOST_DEVICE types::get<axis_reg, L> get_axis() const {
    return get_axis<types::index_cast<axis_reg, L>>();
  }

  /// @tparam axis_t type of the axis.
  /// @returns an axis object of the given type.
  template <typename axis_t>
  DETRAY_HOST_DEVICE axis_t get_axis() const {
    return get_axis<axis_t::bounds_type::label>();
  }
  /// @}

  /// @returns the number of bins per axis
  DETRAY_HOST_DEVICE constexpr auto nbins_per_axis() const -> loc_bin_index {
    // Empty bin indices to be filled
    loc_bin_index n_bins{};
    // Get the number of bins for every axis
    (get_axis_nbins(get_axis<axis_ts>(), n_bins), ...);

    return n_bins;
  }

  /// @returns the total number of bins over all axes
  DETRAY_HOST_DEVICE constexpr auto nbins() const -> dindex {
    const auto n_bins_per_axis = nbins_per_axis();
    dindex n_bins{1u};
    for (dindex i = 0u; i < dim; ++i) {
      n_bins *= n_bins_per_axis[i];
    }
    return n_bins;
  }

  /// Query the bin index for every coordinate of the given point on the axes.
  ///
  /// @param p the point in the local coordinate system that is spanned
  ///          by the axes.
  ///
  /// @returns a multi bin that contains the resulting bin indices for
  ///          every axis in the corresponding entry (e.g. bin_x in entry 0)
  DETRAY_HOST_DEVICE loc_bin_index bins(const point_type &p) const {
    // Empty bin indices to be filled
    loc_bin_index bin_indices{};
    // Run the bin resolution for every axis in this multi-axis type
    (get_axis_bin(get_axis<axis_ts>(), p, bin_indices), ...);

    return bin_indices;
  }

  /// @brief Get a bin range on every axis corresponding to the given point
  /// and the neighborhood around it.
  ///
  /// The neighborhood around the point can be defined in two ways:
  /// - scalar: neighborhood around the axis value,
  /// - index: neighborhood in #bins
  /// The resulting bin index range will contain all bins that belong to a
  /// given neighborhood around the lookup point.
  ///
  /// @tparam neighbor_t the type of neighborhood defined on the axis around
  ///                    the point
  ///
  /// @param p the point in the local coordinate system that is spanned
  ///          by the axes.
  /// @param nhood the search window definition.
  ///
  /// @returns a multi bin range that contains the resulting bin ranges for
  ///          every axis in the corresponding entry (e.g. rng_x in entry 0)
  template <typename neighbor_t>
  DETRAY_HOST_DEVICE multi_bin_range<dim> bin_ranges(
      const point_type &p, const darray<neighbor_t, 2> &nhood) const {
    // Empty bin ranges to be filled
    multi_bin_range<dim> bin_ranges{};
    // Run the range resolution for every axis in this multi-axis type
    (get_axis_bin_ranges(get_axis<axis_ts>(), p, nhood, bin_ranges), ...);

    return bin_ranges;
  }

  /// @returns a vecmem view on the axes data. Only allowed if it owns data.
  template <bool owner = is_owning>
    requires owner
  DETRAY_HOST auto get_data() -> view_type {
    return view_type{detray::get_data(m_edge_offsets),
                     detray::get_data(m_edges)};
  }

  /// @returns a vecmem const view on the axes data. Only allowed if it is
  /// owning data.
  template <bool owner = is_owning>
    requires owner
  DETRAY_HOST auto get_data() const -> const_view_type {
    return const_view_type{detray::get_data(m_edge_offsets),
                           detray::get_data(m_edges)};
  }

  /// Equality operator
  ///
  /// @param rhs the right-hand side of the comparison
  ///
  /// @note in the non-owning case, we compare the values not the pointers
  ///
  /// @returns whether the two axes are equal
  DETRAY_HOST_DEVICE constexpr auto operator==(const multi_axis &rhs) const
      -> bool {
    if constexpr (!std::is_pointer_v<edge_range_t>) {
      return m_edge_offsets == rhs.m_edge_offsets && m_edges == rhs.m_edges;
    } else {
      return m_edge_offsets == rhs.m_edge_offsets && *m_edges == *rhs.m_edges;
    }
    return false;
  }

  /// @returns a string stream that prints the multi axis details
  DETRAY_HOST
  friend std::ostream &operator<<(std::ostream &os, const multi_axis &ax) {
    os << "Axis 0:\n" << ax.template get_axis<0>();

    if constexpr (multi_axis::dim > 1) {
      os << "\nAxis 1:\n" << ax.template get_axis<1>();
    }
    if constexpr (multi_axis::dim > 2) {
      os << "\nAxis 2:\n" << ax.template get_axis<2>();
    }

    return os;
  }

 private:
  /// Get the number of bins for a single axis.
  ///
  /// @tparam axis_t defines the axis for the lookup (axis types are unique)
  ///
  /// @param [in] ax the axis that performs the lookup
  /// @param [out] n_bins the resulting bin numbers
  template <typename axis_t>
  DETRAY_HOST_DEVICE void get_axis_nbins(const axis_t &ax,
                                         loc_bin_index &n_bins) const {
    // Get the index corresponding to the axis label (e.g. bin_x <=> 0)
    constexpr auto loc_idx{
        types::index_cast<axis_reg, axis_t::bounds_type::label>};
    n_bins[loc_idx] = ax.nbins();
  }

  /// Perform the bin lookup on a particular axis
  ///
  /// @tparam axis_t defines the axis for the lookup (axis types are unique)
  ///
  /// @param [in] ax the axis that performs the lookup
  /// @param [in] p the point to be looked up on the axis
  /// @param [out] bin_indices the multi-bin object that is filled with the
  ///                          loc bin indices (axis index corresponds to
  ///                          entry of the multi-bin (e.g. binx <=> 0))
  template <typename axis_t>
  DETRAY_HOST_DEVICE void get_axis_bin(const axis_t &ax, const point_type &p,
                                       loc_bin_index &bin_indices) const {
    // Get the index corresponding to the axis label (e.g. bin_x <=> 0)
    constexpr auto loc_idx{
        types::index_cast<axis_reg, axis_t::bounds_type::label>};
    bin_indices[loc_idx] = ax.bin(p[loc_idx]);
  }

  /// Perform the bin lookup on a particular axis within a given bin
  /// neighborhood
  ///
  /// @tparam axis_t defines the axis for the lookup (axis types are unique)
  /// @tparam neighbor_t the type of neighborhood defined on the axis around
  ///                    the point
  ///
  /// @param [in] ax the axis that performs the lookup
  /// @param [in] p the point to be looked up on the axis
  /// @param [in] nhood the neighborhood around the point for the range lookup
  /// @param [out] bin_ranges the multi-bin-range object that is filled with
  ///                         the neighbor bin range
  template <typename axis_t, typename neighbor_t>
  DETRAY_HOST_DEVICE void get_axis_bin_ranges(
      const axis_t &ax, const point_type &p, const darray<neighbor_t, 2> &nhood,
      multi_bin_range<dim> &bin_ranges) const {
    // Get the index corresponding to the axis label (e.g. bin_range_x = 0)
    constexpr auto loc_idx{
        types::index_cast<axis_reg, axis_t::bounds_type::label>};
    bin_ranges[loc_idx] = ax.range(p[loc_idx], nhood);

    assert(bin_ranges[loc_idx][0] <= bin_ranges[loc_idx][1]);
  }

  /// Data that the axes keep: index ranges in the edges container
  edge_offset_range_t m_edge_offsets{};
  /// Contains all bin edges for all axes
  edge_range_t m_edges{};
};

namespace detail {

/// @brief Helper type to assemble a multi-axis from bounds and binnings
template <bool is_owning, typename containers, typename local_frame, typename,
          typename>
struct multi_axis_assembler;

/// @brief Specialized struct to extract axis bounds and binnings from a tuple
template <bool is_owning, typename containers, typename local_frame,
          typename... axis_bounds, typename... binning_ts>
struct multi_axis_assembler<is_owning, containers, local_frame,
                            types::list<axis_bounds...>,
                            types::list<binning_ts...>> {
  static_assert(sizeof...(axis_bounds) > 0,
                "At least one bounds type needs to be defined");
  static_assert(sizeof...(axis_bounds) == sizeof...(binning_ts),
                "Number of axis bounds for this mask and given binning types "
                "don't match!");

  using type = axis::multi_axis<is_owning, local_frame,
                                axis::single_axis<axis_bounds, binning_ts>...>;
};

}  // namespace detail

}  // namespace detray::axis

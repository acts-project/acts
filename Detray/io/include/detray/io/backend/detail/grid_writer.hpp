// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2023-2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/definitions/indexing.hpp"
#include "detray/io/backend/detail/basic_converter.hpp"
#include "detray/io/backend/detail/type_info.hpp"
#include "detray/io/frontend/payloads.hpp"
#include "detray/utils/detector_statistics.hpp"
#include "detray/utils/grid/concepts.hpp"
#include "detray/utils/grid/grid.hpp"
#include "detray/utils/type_list.hpp"
#include "detray/utils/type_traits.hpp"

// System include(s)
#include <algorithm>
#include <string>
#include <string_view>
#include <type_traits>
#include <vector>

namespace detray::io::detail {

/// @brief Grid writer backend utility
class grid_writer {
 public:
  /// Convert the header information into its payload
  template <typename grid_store_t>
  static grid_header_payload header_to_payload(
      const std::string_view writer_tag, const grid_store_t& store,
      const std::string_view det_name) {
    grid_header_payload header_data;

    header_data.common =
        detail::basic_converter::to_payload(det_name, writer_tag);

    header_data.sub_header.emplace();
    auto& grid_sub_header = header_data.sub_header.value();
    grid_sub_header.n_grids = detray::n_grids(store);

    return header_data;
  }

 protected:
  /// Convert a grid from a collection into its payload
  ///
  /// @param store the data store of grids (tuple of grid collections)
  /// @param grid_link type and index of the grid
  /// @param vol_idx type and index of volume the grid belongs to
  /// @param owner_idx inder of the owner of the grid (can be identical to
  /// vol_idx or volume-local surface index)
  /// @param grid_data the grid payload to be filled
  /// @param converter callable that can convert a grid bin entry into its
  /// respective IO payload (of type @tparam content_t)
  template <typename store_t, typename content_t, typename grid_id_t,
            typename converter_t>
  static void to_payload(
      const store_t& store, typename store_t::single_link grid_link,
      dindex vol_idx, dindex owner_idx,
      detector_grids_payload<content_t, grid_id_t>& grids_data,
      converter_t converter) {
    // If the accelerator is a grid, insert the payload
    store.template visit<get_grid_payload>(grid_link, vol_idx, owner_idx,
                                           grids_data, converter);
  }

  /// Convert a grid @param gr of type @param type and index @param idx
  /// into its io payload, using @param converter for the bin content
  template <typename content_t, typename value_t, typename grid_id_t,
            class grid_t>
  static grid_payload<content_t, grid_id_t> to_payload(
      std::size_t owner_index, grid_id_t type, const std::size_t idx,
      const grid_t& gr, std::function<content_t(const value_t&)> converter) {
    if (type == grid_id_t::unknown) {
      types::print<types::list<grid_t>>();
      throw std::invalid_argument("Could not match type to IO type-id");
    }

    grid_payload<content_t, grid_id_t> grid_data;

    grid_data.owner_link = detail::basic_converter::to_payload(owner_index);
    grid_data.grid_link = detail::basic_converter::to_payload(type, idx);

    // Convert the multi-axis into single axis payloads
    const darray<axis_payload, grid_t::dim> axes_data = to_payload(gr.axes());

    grid_data.axes.resize(axes_data.size());
    std::ranges::copy(axes_data, std::begin(grid_data.axes));

    // Write the surface indices
    for (unsigned int gid = 0u; gid < gr.nbins(); ++gid) {
      // Get the local bin indices and convert the bin into its payload
      grid_bin_payload binp =
          to_payload(gr.deserialize(gid), gr.bin(gid), converter);
      grid_data.bins.push_back(std::move(binp));
    }

    return grid_data;
  }

  /// Convert a multi-axis @param axes into its io payload
  template <bool ownership, typename local_frame_t, typename... axis_ts>
  static auto to_payload(
      const axis::multi_axis<ownership, local_frame_t, axis_ts...>& axes) {
    // Convert every single axis and construct array from their payloads
    darray<axis_payload, sizeof...(axis_ts)> axes_data{
        to_payload(axes.template get_axis<axis_ts>())...};

    return axes_data;
  }

  /// Convert a single axis @param axis into its io payload
  template <typename bounds_t, typename binning_t>
  static axis_payload to_payload(
      const axis::single_axis<bounds_t, binning_t>& axis) {
    axis_payload axis_data;

    axis_data.binning = axis.binning();
    axis_data.bounds = axis.bounds();
    axis_data.label = axis.label();
    axis_data.bins = axis.nbins();

    if (axis.binning() == axis::binning::e_regular) {
      axis_data.edges = {axis.min(), axis.max()};
    } else {
      const auto& bin_edges = axis.bin_edges();
      axis_data.edges.resize(bin_edges.size());
      std::ranges::copy(bin_edges, std::begin(axis_data.edges));
    }

    return axis_data;
  }

  /// Convert a multi-bin @param mbin into its io payload
  template <typename content_t, typename value_t, std::size_t DIM,
            typename content_range_t>
  static grid_bin_payload<content_t> to_payload(
      const axis::multi_bin<DIM> mbin, const content_range_t& content,
      std::function<content_t(const value_t&)> converter) {
    grid_bin_payload<content_t> bin_data;

    // Local bin indices are written in the order the grid axis are stored
    for (unsigned int i = 0u; i < DIM; ++i) {
      bin_data.loc_index.push_back(mbin[i]);
    }

    // Put all entries of the bin into the payload
    bin_data.content.reserve(content.size());
    for (const auto& entry : content) {
      bin_data.content.push_back(converter(entry));
    }

    return bin_data;
  }

 private:
  /// Retrieve a @c grid_payload from grid collection elements
  struct get_grid_payload {
    template <typename grid_group_t, typename index_t, typename content_t,
              typename grid_id_t, typename converter_t>
    inline void operator()(
        [[maybe_unused]] const grid_group_t& coll,
        [[maybe_unused]] const index_t& index,
        [[maybe_unused]] std::size_t vol_idx,
        [[maybe_unused]] std::size_t owner_idx,
        [[maybe_unused]] detector_grids_payload<content_t, grid_id_t>&
            grids_data,
        [[maybe_unused]] converter_t& converter) const {
      using coll_value_t = typename grid_group_t::value_type;
      using value_t = detray::detail::get_value_t<coll_value_t>;

      // Only try to convert grids that match the converter that was
      // passed
      constexpr bool is_convertible{std::invocable<converter_t, value_t>};

      if constexpr (detray::concepts::grid<coll_value_t> && is_convertible) {
        auto gr_pyload = to_payload<content_t, value_t>(
            owner_idx, io::detail::get_id<coll_value_t>(), index, coll[index],
            converter);

        auto& grids_map = grids_data.grids;
        auto search = grids_map.find(vol_idx);
        if (search != grids_map.end()) {
          grids_map.at(vol_idx).push_back(std::move(gr_pyload));
        } else {
          grids_map[vol_idx] = {std::move(gr_pyload)};
        }
      }
    }
  };
};

}  // namespace detray::io::detail

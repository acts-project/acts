// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// Project include(s)
#include "detray/definitions/containers.hpp"
#include "detray/definitions/indexing.hpp"
#include "detray/material/material.hpp"
#include "detray/material/material_rod.hpp"
#include "detray/material/material_slab.hpp"

// System include(s)
#include <algorithm>
#include <cstddef>
#include <functional>
#include <optional>
#include <unordered_map>
#include <utility>

/// Helpers to find identical surface material in the detector material store
/// while building, so that several surfaces can share one material entry.
///
/// @note The comparisons are stricter than the @c operator== of the material
/// types: two entries are only considered identical, if every parameter that
/// enters the material interaction is bitwise equal (up to signed zero), so
/// that sharing an entry cannot change any material lookup.
namespace detray::detail {

/// @returns true if the axes @param lhs and @param rhs have the same type
/// and bin edges
template <typename axis_lhs_t, typename axis_rhs_t>
DETRAY_HOST bool is_identical_axis(const axis_lhs_t &lhs,
                                   const axis_rhs_t &rhs) {
  if (lhs.label() != rhs.label() || lhs.bounds() != rhs.bounds() ||
      lhs.binning() != rhs.binning() || lhs.nbins() != rhs.nbins() ||
      lhs.span() != rhs.span()) {
    return false;
  }
  const auto lhs_edges = lhs.bin_edges();
  const auto rhs_edges = rhs.bin_edges();

  return std::ranges::equal(lhs_edges, rhs_edges);
}

/// @returns true if the material grids @param lhs and @param rhs have
/// identical axes and bin content
///
/// @note Works across owning and non-owning grids. The @c operator== of the
/// grid cannot be used here: For non-owning grids from the same collection it
/// compares the bin edge offsets, i.e. the storage location, not the values
template <typename grid_lhs_t, typename grid_rhs_t>
DETRAY_HOST bool is_identical_grid(const grid_lhs_t &lhs,
                                   const grid_rhs_t &rhs) {
  static_assert(grid_lhs_t::dim == grid_rhs_t::dim);

  if (lhs.nbins() != rhs.nbins()) {
    return false;
  }

  // Compare the axes
  const bool same_axes = [&]<std::size_t... I>(std::index_sequence<I...>) {
    return (is_identical_axis(lhs.template get_axis<I>(),
                              rhs.template get_axis<I>()) &&
            ...);
  }(std::make_index_sequence<grid_lhs_t::dim>{});

  if (!same_axes) {
    return false;
  }

  // Compare the bin content
  for (dindex gbin = 0u; gbin < lhs.nbins(); ++gbin) {
    const auto &lhs_bin = lhs.bin(gbin);
    const auto &rhs_bin = rhs.bin(gbin);

    auto rhs_itr = rhs_bin.begin();
    for (const auto &lhs_entry : lhs_bin) {
      if (rhs_itr == rhs_bin.end() || lhs_entry != *rhs_itr) {
        return false;
      }
      ++rhs_itr;
    }
    if (rhs_itr != rhs_bin.end()) {
      return false;
    }
  }

  return true;
}

}  // namespace detray::detail

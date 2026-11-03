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

/// @returns true if all parameters of the materials @param lhs and @param rhs
/// are equal
template <concepts::scalar scalar_t>
DETRAY_HOST constexpr bool is_identical(const material<scalar_t> &lhs,
                                        const material<scalar_t> &rhs) {
  return lhs == rhs && lhs.mass_density() == rhs.mass_density() &&
         lhs.molar_density() == rhs.molar_density() &&
         lhs.state() == rhs.state() &&
         lhs.has_density_effect_data() == rhs.has_density_effect_data() &&
         lhs.density_effect_data() == rhs.density_effect_data();
}

/// @returns true if the material slabs @param lhs and @param rhs are identical
template <concepts::scalar scalar_t>
DETRAY_HOST constexpr bool is_identical(const material_slab<scalar_t> &lhs,
                                        const material_slab<scalar_t> &rhs) {
  return is_identical(lhs.get_material(), rhs.get_material()) &&
         lhs.thickness() == rhs.thickness() &&
         lhs.thickness_in_X0() == rhs.thickness_in_X0() &&
         lhs.thickness_in_L0() == rhs.thickness_in_L0();
}

/// @returns true if the material rods @param lhs and @param rhs are identical
template <concepts::scalar scalar_t>
DETRAY_HOST constexpr bool is_identical(const material_rod<scalar_t> &lhs,
                                        const material_rod<scalar_t> &rhs) {
  // The radius in X0/L0 is derived from the radius and the material
  return is_identical(lhs.get_material(), rhs.get_material()) &&
         lhs.radius() == rhs.radius();
}

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
      if (rhs_itr == rhs_bin.end() || !is_identical(lhs_entry, *rhs_itr)) {
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

/// @returns a hash of a homogeneous material entry (slab or rod), which is
/// consistent with @c is_identical
template <typename material_entry_t>
DETRAY_HOST std::size_t material_hash(const material_entry_t &entry) {
  using scalar_t = typename material_entry_t::scalar_type;

  const auto &mat = entry.get_material();
  std::size_t seed{0u};
  for (const scalar_t v :
       {entry.thickness(), mat.X0(), mat.L0(), mat.Ar(), mat.Z()}) {
    // Same combination as boost::hash_combine
    seed ^=
        std::hash<scalar_t>{}(v) + 0x9e3779b9u + (seed << 6u) + (seed >> 2u);
  }
  return seed;
}

/// @brief Lookup of identical homogeneous material entries in a collection.
///
/// Holds a hash index over the collection @param coll, which has to outlive
/// the lookup. Entries that are appended to the collection via @c insert are
/// indexed as well.
template <typename collection_t>
class homogeneous_material_lookup {
  using value_type = typename collection_t::value_type;

 public:
  /// Index all entries that are currently in the collection @param coll
  DETRAY_HOST explicit homogeneous_material_lookup(collection_t &coll)
      : m_coll{coll} {
    for (dindex i = 0u; i < static_cast<dindex>(m_coll.size()); ++i) {
      m_index.emplace(material_hash(m_coll[i]), i);
    }
  }

  /// @returns the index of an entry identical to @param entry in the
  /// collection, if one exists
  DETRAY_HOST std::optional<dindex> find(const value_type &entry) const {
    auto [first, last] = m_index.equal_range(material_hash(entry));
    for (auto itr = first; itr != last; ++itr) {
      if (is_identical(m_coll[itr->second], entry)) {
        return itr->second;
      }
    }
    return std::nullopt;
  }

  /// @returns the index of an entry identical to @param entry in the
  /// collection. The entry is appended, if no such entry exists yet.
  DETRAY_HOST dindex insert(const value_type &entry) {
    if (auto idx = find(entry); idx.has_value()) {
      return *idx;
    }
    const auto idx{static_cast<dindex>(m_coll.size())};
    m_coll.push_back(entry);
    m_index.emplace(material_hash(entry), idx);

    return idx;
  }

 private:
  collection_t &m_coll;
  std::unordered_multimap<std::size_t, dindex> m_index{};
};

}  // namespace detray::detail

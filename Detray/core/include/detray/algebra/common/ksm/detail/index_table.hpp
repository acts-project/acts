// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// Project include(s)
#include "detray/algebra/common/ksm/concepts.hpp"
#include "detray/algebra/common/ksm/value.hpp"

// System include(s)
#include <array>
#include <cstddef>
#include <limits>
#include <type_traits>
#include <utility>

// This file assigns storage to the cells of a canonical substructure. Cells
// that name the same index variable share one stored scalar, and a constant
// cell is not stored at all, so the map from a coordinate to a slot in the
// value array is what this file computes.
namespace detray::ksm::detail {
/// @brief Helper functions that get the index of an index variable, or a
/// sentinel value of -1 if another value type is passed.
///
/// @tparam T The type for which we want the index.
///
/// @{
template <typename T>
struct index_variable_getter {
  using type = std::integral_constant<long, -1>;
};

template <std::size_t I>
struct index_variable_getter<index_variable<I>> {
  using type = std::integral_constant<long, I>;
};
/// @}

/// @brief The result of @c make_index_table: the storage slot of every cell,
/// plus the number of distinct slots.
template <std::size_t N>
struct index_table_t {
  std::array<std::size_t, N> table;
  std::size_t num_variables;
};

/// @brief Compute the storage slot of every cell of a canonical substructure.
///
/// Cells are visited in row-major flat order. The first cell to name a given
/// index variable claims the next free slot; a later cell naming the same
/// variable reuses that slot, which is how symmetry and other repeated
/// structure end up sharing one scalar. A cell that is not an index variable
/// is a constant and has no storage, so it is given a sentinel slot.
///
/// A matrix of the following shape (where A..H are non-anonymous variables):
///
/// [[     A,     B,     0],
///  [     D,     E,     F],
///  [     1,     H,     D]]
///
/// Will end up storing a 6-element array in memory, which will look like this:
///
/// [A, B, D, E, F, H]
template <typename Sub, std::size_t... Ks>
  requires(concepts::is_canonical<Sub>)
constexpr auto make_index_table(std::index_sequence<Ks...>) {
  constexpr std::size_t C = Sub::columns;
  constexpr std::array<bool, sizeof...(Ks)> is_indexvar = {
      (concepts::is_index_variable<
          typename Sub::template value_at<Ks / C, Ks % C>>)...};
  constexpr std::array<long, sizeof...(Ks)> indexvar_indices = {
      (index_variable_getter<
          typename Sub::template value_at<Ks / C, Ks % C>>::type::value)...};

  std::array<std::size_t, sizeof...(Ks)> table{};
  std::size_t acc = 0;
  for (std::size_t k = 0; k < sizeof...(Ks); ++k) {
    if (is_indexvar[k]) {
      bool seen_prev = false;
      std::size_t prev_idx = 0;
      for (std::size_t j = 0; j < k; ++j) {
        if (is_indexvar[j] && indexvar_indices[k] == indexvar_indices[j]) {
          seen_prev = true;
          prev_idx = j;
          break;
        }
      }

      if (!seen_prev) {
        table[k] = acc;
        ++acc;
      } else {
        table[k] = table[prev_idx];
      }
    } else {
      table[k] = std::numeric_limits<std::size_t>::max();
    }
  }
  return index_table_t<sizeof...(Ks)>{table, acc};
}

/// @brief The index table of a canonical substructure, and the two projections
/// of it that the substructure itself exposes.
///
/// @{
template <typename Sub>
  requires(concepts::is_canonical<Sub>)
inline constexpr auto index_table_info =
    make_index_table<Sub>(std::make_index_sequence<Sub::rows * Sub::columns>{});

template <typename Sub>
  requires(concepts::is_canonical<Sub>)
inline constexpr auto index_table = index_table_info<Sub>.table;

template <typename Sub>
  requires(concepts::is_canonical<Sub>)
inline constexpr auto get_num_variables = index_table_info<Sub>.num_variables;
/// @}
}  // namespace detray::ksm::detail

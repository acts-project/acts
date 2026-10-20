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
#include "detray/algebra/common/ksm/detail/canonicalization.hpp"
#include "detray/algebra/common/ksm/detail/index_table.hpp"
#include "detray/algebra/common/ksm/detail/operator_substructure.hpp"
#include "detray/algebra/common/ksm/row.hpp"

// System include(s)
#include <cstddef>
#include <tuple>
#include <utility>

// This file defines what a KSM substructure looks like, which contains the
// structural values but contains no storage and no scalar type.
namespace detray::ksm {
namespace detail {
/// @brief Helpers to get a column out of a substructure
///
/// @{
template <typename substructure_t, std::size_t J, typename T>
struct column_getter_helper {};

template <typename substructure_t, std::size_t J, std::size_t... Is>
struct column_getter_helper<substructure_t, J,
                            std::integer_sequence<std::size_t, Is...>> {
  using type = row<typename substructure_t::template value_at<Is, J>...>;
};
/// @}
}  // namespace detail

/// @brief The main substructure class, defined as a list of rows, each of
/// which is itself a list of values.
template <typename... row_ts>
struct substructure {
  using this_t = substructure<row_ts...>;

  // Grab the number of rows, then the number of columns from the first row
  // and assert that all rows are the same size.
  static constexpr std::size_t columns =
      std::tuple_element<0, std::tuple<row_ts...>>::type::num_values;
  static constexpr std::size_t rows = sizeof...(row_ts);
  static_assert(((row_ts::num_values == columns) && ...),
                "Rows do not all have the same size!");

  /// @brief Helper to get a row at a certain index.
  template <std::size_t I>
  using row_at = typename std::tuple_element<I, std::tuple<row_ts...>>::type;

  /// @brief Helper to get the value at a certain index pair.
  template <std::size_t I, std::size_t J>
  using value_at = typename row_at<I>::template value_at<J>;

  /// @brief Compute the canonical version of this substructure, which will
  /// be used for virtually all substructure-related operations.
  using canonical_type =
      typename detail::canonicalize_substructure<this_t>::type;

  /// @brief Get the storage index of a given (I, J) coordinate. Although the
  /// storage of values is not a substructure concern, we keep this here in
  /// order to reuse the computation between matrices with e.g. the same
  /// substructure but different scalar types.
  template <std::size_t I, std::size_t J>
  static constexpr std::size_t element_index =
      detail::index_table<canonical_type>[I * columns + J];

  /// @brief The number of unique variables.
  static constexpr std::size_t num_variables =
      detail::get_num_variables<canonical_type>;

  static_assert((concepts::is_row<row_ts> && ...));

  /// @brief Helper to get the column at a certain index; useful for matrix
  /// multiplication.
  template <std::size_t J>
  using column_at = typename detail::column_getter_helper<
      substructure<row_ts...>, J, std::make_index_sequence<rows>>::type;

  /// @brief Helper to compute the substructure of this + X.
  template <typename other_substructure_t>
    requires(concepts::is_canonical<this_t> &&
             concepts::is_canonical<other_substructure_t>)
  using addition_type =
      detail::additive_substructure<this_t,
                                    other_substructure_t>::type::canonical_type;

  /// @brief Helper to compute the substructure of this - X.
  template <typename other_substructure_t>
    requires(concepts::is_canonical<this_t> &&
             concepts::is_canonical<other_substructure_t>)
  using subtraction_type = detail::subtractive_substructure<
      this_t, other_substructure_t>::type::canonical_type;

  /// @brief Helper to compute the substructure of this * X.
  template <typename other_substructure_t>
    requires(concepts::is_canonical<this_t> &&
             concepts::is_canonical<other_substructure_t>)
  using multiplication_type = detail::multiplicative_substructure<
      this_t, other_substructure_t>::type::canonical_type;

  /// @brief Helper to compute the substructure of this * X * this^T.
  template <typename other_substructure_t>
    requires(concepts::is_canonical<this_t> &&
             concepts::is_canonical<other_substructure_t> &&
             concepts::is_symmetric<other_substructure_t>)
  using congruence_type = detail::congruence_substructure<
      this_t, other_substructure_t>::type::canonical_type;
};
}  // namespace detray::ksm

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
#include "detray/algebra/common/ksm/matrix.hpp"
#include "detray/algebra/common/ksm/value.hpp"
#include "detray/definitions/detail/qualifiers.hpp"

// System include(s).
#include <concepts>
#include <cstddef>
#include <iomanip>
#include <ostream>
#include <sstream>
#include <string>
#include <utility>

namespace detray::ksm {

/// @brief Print a known-substructure matrix
///
/// A cell is addressable at compile time only, so the elements are visited
/// through an index sequence rather than a loop. A structural cell has no
/// storage, but @c at still yields its compile-time value, so the printed
/// matrix is the full one and not only the stored part.
///
/// The kinds of cell are told apart in the output:
///
/// - a plain number is a scalar the matrix stores;
/// - @c &V is a scalar it stores, but shares with an earlier cell, so it is
///   not storage of its own -- the upper triangle of a symmetric matrix is
///   plain and the lower one is marked;
/// - @c . is a structural zero;
/// - @c 'C is any other constant the substructure fixes at compile time.
///
/// A substructure therefore shows both its shape and what it actually costs
/// at a glance: the unmarked cells are exactly the stored scalars.
///
/// @param out the stream to write to
/// @param m the matrix to print
///
/// @returns @p out, so that the calls can be chained
template <typename substructure_t, typename scalar_t>
DETRAY_HOST std::ostream &operator<<(
    std::ostream &out, const matrix<substructure_t, scalar_t> &m) {
  using matrix_t = matrix<substructure_t, scalar_t>;

  constexpr std::size_t n_rows{matrix_t::rows};
  constexpr std::size_t n_cols{matrix_t::columns};

  const auto print_cell = [&out, &m]<std::size_t I, std::size_t J>() {
    // A cell is structural exactly when its canonical value is a constant: an
    // index variable is stored, anything else the substructure fixes.
    using value_t =
        typename matrix_t::canonical_substructure_type::template value_at<I, J>;

    // A stored cell shares its scalar with an earlier one exactly when the
    // substructure names the same index variable in both places, so the cell
    // that first claimed the storage is the one printed unmarked.
    constexpr bool is_shared =
        concepts::is_index_variable<value_t> &&
        detail::first_equal_flat_index_v<
            typename matrix_t::canonical_substructure_type, value_t> !=
            I * n_cols + J;

    const auto emit = [&out, &m](int width) {
      // A marked value is padded as one piece, so that its digits keep the
      // same right edge as an unmarked number in the column above it. How wide
      // the piece is is only known once the scalar is formatted, hence the
      // buffer; @c copyfmt carries over any precision or base the caller set
      // on @p out, so a marked cell formats like an unmarked one.
      const auto marked = [&out, &m](int w, char mark) {
        std::ostringstream buf;
        buf.copyfmt(out);
        buf.width(0);
        buf << mark << m.template at<I, J>();
        out << std::setw(w) << buf.str();
      };

      if constexpr (std::same_as<value_t, zero>) {
        out << std::setw(width) << ".";
      } else if constexpr (concepts::is_integral_value<value_t>) {
        marked(width, '\'');
      } else if constexpr (is_shared) {
        marked(width, '&');
      } else {
        out << std::setw(width) << m.template at<I, J>();
      }
    };

    if constexpr (n_cols == 1u) {
      emit(I == 0u ? 15 : 16);
    } else {
      if constexpr (J == 0u) {
        out << (I == 0u ? "[" : " [");
      }
      emit(J == 0u ? 15 : 16);
      if constexpr (J == n_cols - 1u) {
        out << "]";
        if constexpr (I != n_rows - 1u) {
          out << "\n";
        }
      }
    }
  };

  out << "[";
  [&print_cell]<std::size_t... Ks>(std::index_sequence<Ks...>) {
    (print_cell.template operator()<Ks / n_cols, Ks % n_cols>(), ...);
  }(std::make_index_sequence<n_rows * n_cols>{});
  out << "]";

  return out;
}

}  // namespace detray::ksm

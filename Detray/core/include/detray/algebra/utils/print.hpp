// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// Project include(s)
#include "detray/algebra/concepts.hpp"
#include "detray/definitions/detail/qualifiers.hpp"

// System include(s).
#include <iomanip>
#include <iostream>

namespace detray {

namespace algebra {
/// Print a vector or point as a comma-separated sequence.
/// @param out Stream receiving the formatted value.
/// @param v Vector or point to print.
/// @returns The output stream.
template <typename vector_t>
  requires(concepts::vector<vector_t> || concepts::point<vector_t>)
DETRAY_HOST std::ostream& operator<<(std::ostream& out, const vector_t& v) {
  using index_t = detray::traits::index_t<vector_t>;

  constexpr index_t size{detray::traits::size<vector_t>};

  out << "[";
  for (index_t i = 0; i < size; ++i) {
    out << v[i];
    if (i != size - 1) {
      out << ", ";
    }
  }
  out << "]";

  return out;
}

/// Print a column matrix as a bracketed sequence of rows.
/// @param out Stream receiving the formatted value.
/// @param v Column matrix to print.
/// @returns The output stream.
template <concepts::column_matrix vector_t>
DETRAY_HOST std::ostream& operator<<(std::ostream& out, const vector_t& v) {
  using index_t = detray::traits::index_t<vector_t>;
  using element_getter_t = detray::traits::element_getter_t<vector_t>;

  constexpr index_t rows{detray::traits::rows<vector_t>};

  out << "[";
  for (index_t i = 0; i < rows; ++i) {
    // Account for the sign of negative elements
    out << std::setw(i == 0 ? 15 : 16) << element_getter_t{}(v, i, 0);
  }
  out << "]";

  return out;
}

/// Print a matrix as rows of bracketed elements.
/// @param out Stream receiving the formatted value.
/// @param m Matrix to print.
/// @returns The output stream.
template <concepts::matrix matrix_t>
DETRAY_HOST std::ostream& operator<<(std::ostream& out, const matrix_t& m) {
  using index_t = detray::traits::index_t<matrix_t>;
  using element_getter_t = detray::traits::element_getter_t<matrix_t>;

  constexpr index_t rows{detray::traits::rows<matrix_t>};
  constexpr index_t columns{detray::traits::columns<matrix_t>};

  out << "[";
  for (index_t i = 0; i < rows; ++i) {
    out << (i == 0 ? "[" : " [");
    for (index_t j = 0; j < columns; ++j) {
      // Account for the sign of negative elements
      out << std::setw(j == 0 ? 15 : 16) << element_getter_t{}(m, i, j);
    }
    out << "]";
    if (i != rows - 1) {
      out << "\n";
    }
  }
  out << "]";

  return out;
}

/// Print the matrix of a three-dimensional transform.
/// @param out Stream receiving the formatted value.
/// @param trf Transform to print.
/// @returns The output stream.
template <concepts::transform3D transform_t>
DETRAY_HOST std::ostream& operator<<(std::ostream& out,
                                     const transform_t& trf) {
  out << trf.matrix();

  return out;
}

}  // namespace algebra

// Pull in the print operator definitions for the algebra types
using algebra::operator<<;

}  // namespace detray

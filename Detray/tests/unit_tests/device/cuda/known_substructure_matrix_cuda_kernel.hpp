// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// Detray core include(s)
#include "detray/definitions/detail/qualifiers.hpp"
#include "detray/utils/known_substructure_matrix.hpp"

// System include(s)
#include <cstddef>
#include <utility>

namespace detray::ksm_cuda_test {

// `detray::matrix` is a namespace, so `ksm::matrix` has to stay qualified in
// here -- which is what a caller inside `namespace detray` has to do anyway.
// The cell types collide with nothing and are pulled in for readability.
using ksm::one;
using ksm::row;
using ksm::substructure;
using ksm::variable;
using ksm::zero;

/// The substructure of the full Jacobian, taken from the assertions that
/// transport_covariance_to_bound_impl already emits: 25 free values, two
/// structural ones, nine structural zeros.
using j_full = substructure<
    row<variable, variable, variable, variable, variable, zero>,
    row<variable, variable, variable, variable, variable, zero>,
    row<variable, variable, variable, variable, variable, zero>,
    row<variable, variable, variable, variable, variable, zero>,
    row<zero, zero, zero, zero, one, zero>,
    row<variable, variable, variable, variable, variable, one>>::canonical_type;

/// A symmetric 6x6, standing in for a track covariance. Needed to reach
/// @c congruence, whose result type is declared symmetric.
using sym6 = ksm::make_symmetric_substructure<6>::canonical_type;

/// A dense 2x2 -- the largest shape @c inverse() supports. Inverting it goes
/// through @c inverted_substructure, whose accessor uses a lambda in an
/// unevaluated context: the construct in the header most likely to trip nvcc,
/// and therefore the one most worth exercising here.
using dense2 = ksm::make_dense_substructure<2, 2>::canonical_type;

/// Number of scalars @c compute writes: seven 6x6 results and one 2x2.
inline constexpr std::size_t result_size = 7u * 6u * 6u + 2u * 2u;

/// Write every cell of @p m -- structural cells included -- to @p out.
/// @{
template <typename matrix_t, std::size_t... Ks>
DETRAY_HOST_DEVICE inline void emit(const matrix_t &m, double *out,
                                    std::index_sequence<Ks...>) {
  ((out[Ks] = m.template at<Ks / matrix_t::columns, Ks % matrix_t::columns>()),
   ...);
}

/// @returns the write cursor advanced past the cells of @p m
template <typename matrix_t>
DETRAY_HOST_DEVICE inline double *emit(const matrix_t &m, double *out) {
  emit(m, out, std::make_index_sequence<matrix_t::rows * matrix_t::columns>{});
  return out + matrix_t::rows * matrix_t::columns;
}
/// @}

/// Give every stored cell of @p m a distinct, deterministic value. Structural
/// cells have no storage and are skipped.
/// @{
template <typename matrix_t, std::size_t I, std::size_t J>
DETRAY_HOST_DEVICE inline void set_cell(matrix_t &m, double base) {
  using value_t =
      typename matrix_t::canonical_substructure_type::template value_at<I, J>;

  if constexpr (ksm::detail::checks::is_index_variable<value_t>) {
    m.template at<I, J>() = base + static_cast<double>(value_t::index);
  }
}

template <typename matrix_t, std::size_t... Ks>
DETRAY_HOST_DEVICE inline void fill(matrix_t &m, double base,
                                    std::index_sequence<Ks...>) {
  (set_cell<matrix_t, Ks / matrix_t::columns, Ks % matrix_t::columns>(m, base),
   ...);
}

template <typename matrix_t>
DETRAY_HOST_DEVICE inline void fill(matrix_t &m, double base) {
  fill(m, base, std::make_index_sequence<matrix_t::rows * matrix_t::columns>{});
}
/// @}

/// Exercise every runtime member of @c ksm::matrix and write the results cell
/// by cell into @p out, which must have room for @c result_size scalars.
///
/// Host and device both compile this one definition, so the test can compare
/// the two outputs and any divergence is a genuine device-code problem rather
/// than a difference in what was measured.
DETRAY_HOST_DEVICE inline void compute(double *out) {
  ksm::matrix<j_full, double> a;
  ksm::matrix<j_full, double> b;
  ksm::matrix<sym6, double> c;
  ksm::matrix<dense2, double> d;

  fill(a, 1.);
  fill(b, 1000.);
  fill(c, 7.);
  // Cells 3, 4, 5, 6 -> determinant 3*6 - 5*4 = -2, so d is invertible.
  fill(d, 3.);

  out = emit(a, out);                // at<I, J>() const and the mutable at()
  out = emit(a * b, out);            // operator* -> fill_with, product_dot
  out = emit(a + b, out);            // operator+
  out = emit(a - b, out);            // binary operator-
  out = emit(-a, out);               // unary operator-
  out = emit(a.transpose(), out);    // transpose
  out = emit(a.congruence(c), out);  // congruence
  out = emit(d.inverse(), out);      // inverse -> inverted_substructure
}

}  // namespace detray::ksm_cuda_test

namespace detray {

/// Run @c ksm_cuda_test::compute in a CUDA kernel, writing
/// @c ksm_cuda_test::result_size scalars into @p out, which must be
/// device-accessible.
void ksm_run_on_device(double *out);

}  // namespace detray

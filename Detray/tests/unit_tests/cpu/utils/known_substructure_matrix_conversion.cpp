// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// Project include(s).
#include "detray/utils/known_substructure_matrix.hpp"

// Test include(s).
#include "detray/test/cpu/algebra_fixture.hpp"

// GoogleTest include(s).
#include <gtest/gtest.h>

// System include(s).
#include <cstddef>
#include <type_traits>
#include <utility>

/// @file
/// Conversions between @c ksm::matrix and a dense detray matrix. Unlike the
/// main known-substructure test, which is algebra agnostic and gated behind
/// DETRAY_BUILD_KNOWN_SUBSTRUCTURE_MATRIX_TEST because it is expensive to
/// compile, these run for every algebra plugin: @c to_dense and @c from_dense
/// are templated on the algebra, so the backends are what they most need
/// covering against.

namespace detray::test {

namespace {

using ksm::one;
using ksm::row;
using ksm::substructure;
using ksm::variable;
using ksm::zero;

/// The substructure of the full Jacobian, transcribed from the assertions that
/// @c transport_covariance_to_bound_impl already emits: 25 free values, two
/// structural ones and nine structural zeros.
using j_full = substructure<
    row<variable, variable, variable, variable, variable, zero>,
    row<variable, variable, variable, variable, variable, zero>,
    row<variable, variable, variable, variable, variable, zero>,
    row<variable, variable, variable, variable, variable, zero>,
    row<zero, zero, zero, zero, one, zero>,
    row<variable, variable, variable, variable, variable, one>>::canonical_type;

/// Mirrored cells name the same variable, so they share one stored value. This
/// is what exercises the shared-cell branch of @c from_dense.
using sym6 = ksm::make_symmetric_substructure<6>::canonical_type;

/// Apply @p f to every cell index of @p Sub at compile time.
template <typename Sub, typename F>
void for_each_cell(F f) {
  [&f]<std::size_t... Ks>(std::index_sequence<Ks...>) {
    (f.template operator()<Ks / Sub::columns, Ks % Sub::columns>(), ...);
  }(std::make_index_sequence<Sub::rows * Sub::columns>{});
}

/// @returns a matrix whose free cells hold small distinct integers.
///
/// Integers keep every comparison below exact: the products stay well inside
/// the range where a float represents integers exactly, so a structure-aware
/// result and a dense one that differ only in which zero terms they skipped
/// must agree bit for bit. That makes @c EXPECT_EQ the right assertion and
/// avoids picking a tolerance.
template <typename Sub, typename scalar_t>
ksm::matrix<Sub, scalar_t> filled(scalar_t base) {
  ksm::matrix<Sub, scalar_t> m;

  for_each_cell<Sub>([&m, base]<std::size_t I, std::size_t J>() {
    using value_t = typename Sub::template value_at<I, J>;

    if constexpr (ksm::detail::checks::is_index_variable<value_t>) {
      m.template at<I, J>() = base + static_cast<scalar_t>(value_t::index);
    }
  });

  return m;
}

}  // namespace

/// Structural cells have no storage, so to_dense has to reconstruct them from
/// their compile-time values.
TEST_F(detray_algebra, ksm_to_dense_materialises_structure) {
  const auto m = filled<j_full>(scalar{1});
  const auto d = m.to_dense<algebra_t>();

  static_assert(
      std::is_same_v<std::decay_t<decltype(d)>, dmatrix<algebra_t, 6, 6>>);

  // The two structural ones ...
  EXPECT_EQ((getter::element<4, 4>(d)), scalar{1});
  EXPECT_EQ((getter::element<5, 5>(d)), scalar{1});
  // ... and a sample of the nine structural zeros.
  EXPECT_EQ((getter::element<0, 5>(d)), scalar{0});
  EXPECT_EQ((getter::element<4, 0>(d)), scalar{0});

  // Every cell, structural or not, must read back what the ksm matrix says.
  for_each_cell<j_full>([&d, &m]<std::size_t I, std::size_t J>() {
    EXPECT_EQ((getter::element<I, J>(d)), (m.template at<I, J>()));
  });
}

/// ksm -> dense -> ksm is the identity.
TEST_F(detray_algebra, ksm_round_trip) {
  const auto m = filled<j_full>(scalar{1});
  const auto back = ksm::matrix<j_full, scalar>::from_dense<algebra_t>(
      m.to_dense<algebra_t>());

  for_each_cell<j_full>([&back, &m]<std::size_t I, std::size_t J>() {
    EXPECT_EQ((back.template at<I, J>()), (m.template at<I, J>()));
  });
}

/// The same, for a substructure whose mirrored cells share storage: from_dense
/// writes each shared value once and checks that the rest agree.
TEST_F(detray_algebra, ksm_round_trip_shared_cells) {
  static_assert(sym6::num_variables == 21u,
                "a symmetric 6x6 stores its upper triangle");

  const auto m = filled<sym6>(scalar{1});
  const auto back =
      ksm::matrix<sym6, scalar>::from_dense<algebra_t>(m.to_dense<algebra_t>());

  for_each_cell<sym6>([&back, &m]<std::size_t I, std::size_t J>() {
    EXPECT_EQ((back.template at<I, J>()), (m.template at<I, J>()));
  });
}

/// The payoff check: skipping the structurally zero terms must not change the
/// answer. The structure-aware product has to equal the dense one computed by
/// the algebra plugin.
TEST_F(detray_algebra, ksm_product_matches_dense) {
  const auto a = filled<j_full>(scalar{1});
  const auto b = filled<j_full>(scalar{100});

  const auto structured = (a * b).to_dense<algebra_t>();
  const dmatrix<algebra_t, 6, 6> dense =
      a.to_dense<algebra_t>() * b.to_dense<algebra_t>();

  for_each_cell<j_full>([&structured, &dense]<std::size_t I, std::size_t J>() {
    EXPECT_EQ((getter::element<I, J>(structured)),
              (getter::element<I, J>(dense)));
  });
}

/// Congruence computes only the upper triangle of a result it knows to be
/// symmetric, so it has the most to get wrong. It must still equal the plain
/// triple product.
TEST_F(detray_algebra, ksm_congruence_matches_dense) {
  const auto a = filled<j_full>(scalar{1});
  const auto c = filled<sym6>(scalar{1});

  const auto structured = a.congruence(c).to_dense<algebra_t>();

  const auto a_dense = a.to_dense<algebra_t>();
  const dmatrix<algebra_t, 6, 6> dense =
      a_dense * c.to_dense<algebra_t>() * detray::matrix::transpose(a_dense);

  for_each_cell<j_full>([&structured, &dense]<std::size_t I, std::size_t J>() {
    EXPECT_EQ((getter::element<I, J>(structured)),
              (getter::element<I, J>(dense)));
  });
}

}  // namespace detray::test

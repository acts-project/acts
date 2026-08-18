// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// Project include(s).
#include "detray/propagator/detail/full_jacobian.hpp"

#include "detray/utils/known_substructure_matrix.hpp"
#include "detray/utils/known_substructure_matrix_types.hpp"

// Test include(s).
#include "detray/test/cpu/algebra_fixture.hpp"

// GoogleTest include(s).
#include <gtest/gtest.h>

// System include(s).
#include <cstddef>
#include <utility>

/// @file
/// The Jacobian products that used to be sympy-generated are now carried out
/// through the matrices' substructures. What the generators guaranteed by
/// construction -- that the elided terms really were zero -- has to be
/// guaranteed by a test instead, so each product here is checked against the
/// same product evaluated densely by the algebra plugin.
///
/// The fills are small integers on purpose. Every intermediate then stays well
/// inside the range where a float represents integers exactly, so a
/// structure-aware result and a dense one that differ only in which zero terms
/// they skipped must agree bit for bit. That makes EXPECT_EQ the right
/// assertion, and removes the need to justify a tolerance.

namespace detray::test {

namespace {

/// Apply @p f to every cell index of an R x C matrix at compile time.
template <std::size_t R, std::size_t C, typename F>
void for_each_cell(F f) {
  [&f]<std::size_t... Ks>(std::index_sequence<Ks...>) {
    (f.template operator()<Ks / C, Ks % C>(), ...);
  }(std::make_index_sequence<R * C>{});
}

/// Fill a dense matrix with small distinct integers.
template <typename matrix_t, typename scalar_t>
void fill(matrix_t &m, std::size_t rows, std::size_t cols, scalar_t base) {
  for (std::size_t i = 0u; i < rows; ++i) {
    for (std::size_t j = 0u; j < cols; ++j) {
      getter::element(m, i, j) =
          base + static_cast<scalar_t>((i * cols + j) % 5u);
    }
  }
}

}  // namespace

/// The full Jacobian: J_F2B * (D + I) * J_transport * J_B2F, where the factors'
/// zeros are elided, against the same product taken densely.
///
/// All four frame/gradient combinations are covered. The two booleans are the
/// frame-dependent zeros: only line frames give d(pos)/d(angle) a value or set
/// the direction terms of the path derivative.
template <bool has_dpos_dangle, bool path_has_direction_terms, bool gradient,
          typename algebra_t, typename scalar_t>
void check_full_jacobian() {
  using transport_t =
      ksm::matrix<ksm::transport_jacobian_substructure<gradient>, scalar_t>;

  dmatrix<algebra_t, 3, 2> b2f_dpos_dloc{}, b2f_ddir_dangle{},
      b2f_dpos_dangle{};
  dmatrix<algebra_t, 2, 3> f2b_dloc_dpos{}, f2b_dangle_ddir{};
  dmatrix<algebra_t, 8, 1> p2f{};
  dmatrix<algebra_t, 1, 8> f2p{};

  fill(b2f_dpos_dloc, 3u, 2u, scalar_t{1});
  fill(b2f_ddir_dangle, 3u, 2u, scalar_t{2});
  fill(f2b_dloc_dpos, 2u, 3u, scalar_t{1});
  fill(f2b_dangle_ddir, 2u, 3u, scalar_t{2});
  fill(p2f, 8u, 1u, scalar_t{1});
  fill(f2p, 1u, 8u, scalar_t{1});

  // The structural zeros the substructures declare. update_full_jacobian
  // asserts these hold, so getting them wrong here is itself a check.
  getter::element(b2f_ddir_dangle, 2u, 0u) = scalar_t{0};  // d(n_z)/d(phi)
  getter::element(f2b_dangle_ddir, 0u, 2u) = scalar_t{0};  // d(phi)/d(n_z)
  getter::element(p2f, 3u, 0u) = scalar_t{0};              // free time
  getter::element(f2p, 0u, 3u) = scalar_t{0};
  getter::element(f2p, 0u, 7u) = scalar_t{0};
  if constexpr (has_dpos_dangle) {
    fill(b2f_dpos_dangle, 3u, 2u, scalar_t{3});
  }
  if constexpr (!path_has_direction_terms) {
    for (std::size_t i = 4u; i < 7u; ++i) {
      getter::element(f2p, 0u, i) = scalar_t{0};
    }
  }

  transport_t transport = transport_t::identity();
  for_each_cell<8u, 8u>([&transport]<std::size_t I, std::size_t J>() {
    using value_t =
        typename transport_t::canonical_substructure_type::template value_at<I,
                                                                             J>;
    if constexpr (ksm::detail::checks::is_index_variable<value_t>) {
      transport.template at<I, J>() =
          scalar_t{1} + static_cast<scalar_t>(value_t::index % 4u);
    }
  });

  // --- the structure-aware product
  const auto structured =
      detail::update_full_jacobian<has_dpos_dangle, path_has_direction_terms,
                                   algebra_t>(
          transport, b2f_dpos_dloc, b2f_ddir_dangle, b2f_dpos_dangle, p2f, f2p,
          f2b_dloc_dpos, f2b_dangle_ddir)
          .template to_dense<algebra_t>();

  // --- the same product, densely
  dmatrix<algebra_t, 8, 6> b2f_dense{};
  dmatrix<algebra_t, 6, 8> f2b_dense{};
  for (std::size_t i = 0u; i < 3u; ++i) {
    for (std::size_t j = 0u; j < 2u; ++j) {
      getter::element(b2f_dense, i, j) = getter::element(b2f_dpos_dloc, i, j);
      getter::element(b2f_dense, i, j + 2u) =
          getter::element(b2f_dpos_dangle, i, j);
      getter::element(b2f_dense, i + 4u, j + 2u) =
          getter::element(b2f_ddir_dangle, i, j);
    }
  }
  getter::element(b2f_dense, 7u, 4u) = scalar_t{1};
  getter::element(b2f_dense, 3u, 5u) = scalar_t{1};
  for (std::size_t i = 0u; i < 2u; ++i) {
    for (std::size_t j = 0u; j < 3u; ++j) {
      getter::element(f2b_dense, i, j) = getter::element(f2b_dloc_dpos, i, j);
      getter::element(f2b_dense, i + 2u, j + 4u) =
          getter::element(f2b_dangle_ddir, i, j);
    }
  }
  getter::element(f2b_dense, 4u, 7u) = scalar_t{1};
  getter::element(f2b_dense, 5u, 3u) = scalar_t{1};

  const dmatrix<algebra_t, 8, 8> d_plus_i =
      p2f * f2p + detray::matrix::identity<dmatrix<algebra_t, 8, 8>>();
  const dmatrix<algebra_t, 6, 6> dense =
      f2b_dense * d_plus_i * transport.template to_dense<algebra_t>() *
      b2f_dense;

  for_each_cell<6u, 6u>([&structured, &dense]<std::size_t I, std::size_t J>() {
    EXPECT_EQ((getter::element<I, J>(structured)),
              (getter::element<I, J>(dense)));
  });
}

TEST_F(detray_algebra, full_jacobian_planar_frames_no_gradient) {
  check_full_jacobian<false, false, false, algebra_t, scalar>();
}

TEST_F(detray_algebra, full_jacobian_planar_frames_with_gradient) {
  check_full_jacobian<false, false, true, algebra_t, scalar>();
}

TEST_F(detray_algebra, full_jacobian_line_frames_no_gradient) {
  check_full_jacobian<true, true, false, algebra_t, scalar>();
}

TEST_F(detray_algebra, full_jacobian_line_frames_with_gradient) {
  check_full_jacobian<true, true, true, algebra_t, scalar>();
}

/// The RK stepper's transport Jacobian step, D * J, likewise carried out in the
/// substructure. Checked against the dense product of the same two matrices.
template <bool gradient, typename algebra_t, typename scalar_t>
void check_transport_step() {
  using transport_t =
      ksm::matrix<ksm::transport_jacobian_substructure<gradient>, scalar_t>;

  transport_t d{}, j{};
  // give both matrices' free cells distinct values
  for_each_cell<8u, 8u>([&d, &j]<std::size_t I, std::size_t J>() {
    using value_t =
        typename transport_t::canonical_substructure_type::template value_at<I,
                                                                             J>;
    if constexpr (ksm::detail::checks::is_index_variable<value_t>) {
      d.template at<I, J>() =
          scalar_t{1} + static_cast<scalar_t>(value_t::index % 4u);
      j.template at<I, J>() =
          scalar_t{2} + static_cast<scalar_t>(value_t::index % 3u);
    }
  });

  const auto structured = (d * j).template to_dense<algebra_t>();
  const dmatrix<algebra_t, 8, 8> dense =
      d.template to_dense<algebra_t>() * j.template to_dense<algebra_t>();

  for_each_cell<8u, 8u>([&structured, &dense]<std::size_t I, std::size_t J>() {
    EXPECT_EQ((getter::element<I, J>(structured)),
              (getter::element<I, J>(dense)));
  });
}

TEST_F(detray_algebra, transport_jacobian_step_no_gradient) {
  check_transport_step<false, algebra_t, scalar>();
}

TEST_F(detray_algebra, transport_jacobian_step_with_gradient) {
  check_transport_step<true, algebra_t, scalar>();
}

/// copy_elements_into has to write the free cells and leave the structural ones
/// to the substructure. Its guard is that a source disagreeing with a
/// structural cell is caught rather than silently dropped.
TEST_F(detray_algebra, copy_elements_into_respects_structure) {
  using b2f_t = ksm::matrix<ksm::bound_to_free_substructure<false>, scalar>;

  dmatrix<algebra_t, 3, 2> ddir_dangle{};
  fill(ddir_dangle, 3u, 2u, scalar{1});
  // (2, 0) lands on a structural zero of the destination
  getter::element(ddir_dangle, 2u, 0u) = scalar{0};

  b2f_t b2f{};
  ksm::copy_elements_into<4u, 2u>(b2f, ddir_dangle);

  for_each_cell<3u, 2u>([&b2f, &ddir_dangle]<std::size_t I, std::size_t J>() {
    EXPECT_EQ((b2f.template at<4u + I, 2u + J>()),
              (getter::element<I, J>(ddir_dangle)));
  });

  // The structural cells the copy did not cover keep their fixed values.
  EXPECT_EQ((b2f.template at<7u, 4u>()), scalar{1});
  EXPECT_EQ((b2f.template at<3u, 5u>()), scalar{1});
  EXPECT_EQ((b2f.template at<0u, 5u>()), scalar{0});
}

/// identity() must reconstruct the identity through a substructure that cannot
/// store most of it.
TEST_F(detray_algebra, ksm_identity) {
  using transport_t =
      ksm::matrix<ksm::transport_jacobian_substructure<false>, scalar>;

  const auto id = transport_t::identity().template to_dense<algebra_t>();

  for_each_cell<8u, 8u>([&id]<std::size_t I, std::size_t J>() {
    EXPECT_EQ((getter::element<I, J>(id)), (I == J) ? scalar{1} : scalar{0});
  });
}

}  // namespace detray::test

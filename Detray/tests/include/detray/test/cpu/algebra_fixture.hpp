// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// Project include(s).
#include "detray/algebra/concepts.hpp"
#include "detray/algebra/type_traits.hpp"

// Test include(s).
#include "detray/test/framework/fixture_base.hpp"
#include "detray/test/framework/types.hpp"

namespace detray::test {

/// Test case class, contains the basic test definitions and algebra types
class detray_algebra : public fixture_base<> {
  using base = fixture_base<>;

 public:
  /// Constructor
  using base::base;

 protected:
#if !DETRAY_ALGEBRA_VC_AOS
  template <std::size_t ROWS, std::size_t COLS>
  void matrix_test_ops_any_matrix() {
    // Test the set_product method.
    {
      matrix<ROWS, ROWS> m1;
      matrix<ROWS, COLS> m2;

      for (std::size_t i = 0; i < ROWS; ++i) {
        for (std::size_t j = 0; j < ROWS; ++j) {
          detray::getter::element(m1, i, j) = static_cast<scalar>(i * ROWS + j);
        }
      }

      for (std::size_t i = 0; i < ROWS; ++i) {
        for (std::size_t j = 0; j < COLS; ++j) {
          detray::getter::element(m2, i, j) = static_cast<scalar>(i * COLS + j);
        }
      }

      {
        matrix<ROWS, COLS> r1 = m1 * m2;
        matrix<ROWS, COLS> r2;
        detray::matrix::set_product(r2, m1, m2);

        for (std::size_t i = 0; i < ROWS; ++i) {
          for (std::size_t j = 0; j < COLS; ++j) {
            ASSERT_NEAR(detray::getter::element(r1, i, j),
                        detray::getter::element(r2, i, j), this->epsilon());
          }
        }
      }
    }

    // Test the set_product_right_transpose method.
    {
      matrix<ROWS, ROWS> m1;
      matrix<COLS, ROWS> m2;

      for (std::size_t i = 0; i < ROWS; ++i) {
        for (std::size_t j = 0; j < ROWS; ++j) {
          detray::getter::element(m1, i, j) = static_cast<scalar>(i * ROWS + j);
        }
      }

      for (std::size_t i = 0; i < COLS; ++i) {
        for (std::size_t j = 0; j < ROWS; ++j) {
          detray::getter::element(m2, i, j) = static_cast<scalar>(i * COLS + j);
        }
      }

      {
        matrix<ROWS, COLS> r1 = m1 * detray::matrix::transpose(m2);
        matrix<ROWS, COLS> r2;
        detray::matrix::set_product_right_transpose(r2, m1, m2);

        for (std::size_t i = 0; i < ROWS; ++i) {
          for (std::size_t j = 0; j < COLS; ++j) {
            ASSERT_NEAR(detray::getter::element(r1, i, j),
                        detray::getter::element(r2, i, j), this->epsilon());
          }
        }
      }
    }

    // Test the set_product_left_transpose method.
    {
      matrix<ROWS, ROWS> m1;
      matrix<ROWS, COLS> m2;

      for (std::size_t i = 0; i < ROWS; ++i) {
        for (std::size_t j = 0; j < ROWS; ++j) {
          detray::getter::element(m1, i, j) = static_cast<scalar>(i * ROWS + j);
        }
      }

      for (std::size_t i = 0; i < ROWS; ++i) {
        for (std::size_t j = 0; j < COLS; ++j) {
          detray::getter::element(m2, i, j) = static_cast<scalar>(i * COLS + j);
        }
      }

      {
        matrix<ROWS, COLS> r1 = detray::matrix::transpose(m1) * m2;
        matrix<ROWS, COLS> r2;
        detray::matrix::set_product_left_transpose(r2, m1, m2);

        for (std::size_t i = 0; i < ROWS; ++i) {
          for (std::size_t j = 0; j < COLS; ++j) {
            ASSERT_NEAR(detray::getter::element(r1, i, j),
                        detray::getter::element(r2, i, j), this->epsilon());
          }
        }
      }

      // Test the transposable_product method.
      {
        // Both untransposed
        {
          matrix<ROWS, COLS> r1 = m1 * m2;
          matrix<ROWS, COLS> r2 =
              detray::matrix::transposed_product<false, false>(m1, m2);

          for (std::size_t i = 0; i < ROWS; ++i) {
            for (std::size_t j = 0; j < COLS; ++j) {
              ASSERT_NEAR(detray::getter::element(r1, i, j),
                          detray::getter::element(r2, i, j), this->epsilon());
            }
          }
        }
      }
    }
  }

  template <std::size_t N>
  void matrix_test_ops_square_matrix() {
    {
      matrix<N, N> m1;
      matrix<N, N> m2;

      for (std::size_t i = 0; i < N; ++i) {
        for (std::size_t j = 0; j < N; ++j) {
          detray::getter::element(m1, i, j) = static_cast<scalar>(i * N + j);
          detray::getter::element(m2, i, j) =
              -1.f * static_cast<scalar>(i * N + j) + 42;
        }
      }

      // Test the set_product method.
      {
        matrix<N, N> r1 = m1 * m2;
        matrix<N, N> r2;
        detray::matrix::set_product(r2, m1, m2);

        for (std::size_t i = 0; i < N; ++i) {
          for (std::size_t j = 0; j < N; ++j) {
            ASSERT_NEAR(detray::getter::element(r1, i, j),
                        detray::getter::element(r2, i, j), this->epsilon());
          }
        }
      }

      // Test the set_product_right_transpose method.
      {
        matrix<N, N> r1 = m1 * detray::matrix::transpose(m2);
        matrix<N, N> r2;
        detray::matrix::set_product_right_transpose(r2, m1, m2);

        for (std::size_t i = 0; i < N; ++i) {
          for (std::size_t j = 0; j < N; ++j) {
            ASSERT_NEAR(detray::getter::element(r1, i, j),
                        detray::getter::element(r2, i, j), this->epsilon());
          }
        }
      }

      // Test the set_product_left_transpose method.
      {
        matrix<N, N> r1 = detray::matrix::transpose(m1) * m2;
        matrix<N, N> r2;
        detray::matrix::set_product_left_transpose(r2, m1, m2);

        for (std::size_t i = 0; i < N; ++i) {
          for (std::size_t j = 0; j < N; ++j) {
            ASSERT_NEAR(detray::getter::element(r1, i, j),
                        detray::getter::element(r2, i, j), this->epsilon());
          }
        }
      }

      // Test the set_inplace_product_right method.
      {
        matrix<N, N> r1 = m1 * m2;
        matrix<N, N> r2 = m1;
        detray::matrix::set_inplace_product_right(r2, m2);

        for (std::size_t i = 0; i < N; ++i) {
          for (std::size_t j = 0; j < N; ++j) {
            ASSERT_NEAR(detray::getter::element(r1, i, j),
                        detray::getter::element(r2, i, j), this->epsilon());
          }
        }
      }

      // Test the set_inplace_product_left method.
      {
        matrix<N, N> r1 = m1 * m2;
        matrix<N, N> r2 = m2;
        detray::matrix::set_inplace_product_left(r2, m1);

        for (std::size_t i = 0; i < N; ++i) {
          for (std::size_t j = 0; j < N; ++j) {
            ASSERT_NEAR(detray::getter::element(r1, i, j),
                        detray::getter::element(r2, i, j), this->epsilon());
          }
        }
      }

      // Test the set_inplace_product_right_transpose method.
      {
        matrix<N, N> r1 = m1 * detray::matrix::transpose(m2);
        matrix<N, N> r2 = m1;
        detray::matrix::set_inplace_product_right_transpose(r2, m2);

        for (std::size_t i = 0; i < N; ++i) {
          for (std::size_t j = 0; j < N; ++j) {
            ASSERT_NEAR(detray::getter::element(r1, i, j),
                        detray::getter::element(r2, i, j), this->epsilon());
          }
        }
      }

      // Test the set_inplace_product_left_transpose method.
      {
        matrix<N, N> r1 = detray::matrix::transpose(m1) * m2;
        matrix<N, N> r2 = m2;
        detray::matrix::set_inplace_product_left_transpose(r2, m1);

        for (std::size_t i = 0; i < N; ++i) {
          for (std::size_t j = 0; j < N; ++j) {
            ASSERT_NEAR(detray::getter::element(r1, i, j),
                        detray::getter::element(r2, i, j), this->epsilon());
          }
        }
      }

      // Test the transposable_product method.
      {
        // Only left transposed
        {
          matrix<N, N> r1 = detray::matrix::transpose(m1) * m2;
          matrix<N, N> r2 =
              detray::matrix::transposed_product<true, false>(m1, m2);

          for (std::size_t i = 0; i < N; ++i) {
            for (std::size_t j = 0; j < N; ++j) {
              ASSERT_NEAR(detray::getter::element(r1, i, j),
                          detray::getter::element(r2, i, j), this->epsilon());
            }
          }
        }

        // Only right transposed
        {
          matrix<N, N> r1 = m1 * detray::matrix::transpose(m2);
          matrix<N, N> r2 =
              detray::matrix::transposed_product<false, true>(m1, m2);

          for (std::size_t i = 0; i < N; ++i) {
            for (std::size_t j = 0; j < N; ++j) {
              ASSERT_NEAR(detray::getter::element(r1, i, j),
                          detray::getter::element(r2, i, j), this->epsilon());
            }
          }
        }

        // Both transposed
        {
          matrix<N, N> r1 =
              detray::matrix::transpose(m1) * detray::matrix::transpose(m2);
          matrix<N, N> r2 =
              detray::matrix::transposed_product<true, true>(m1, m2);

          for (std::size_t i = 0; i < N; ++i) {
            for (std::size_t j = 0; j < N; ++j) {
              ASSERT_NEAR(detray::getter::element(r1, i, j),
                          detray::getter::element(r2, i, j), this->epsilon());
            }
          }
        }
      }
    }

    this->template matrix_test_ops_any_matrix<N, N>();
  }

  template <std::size_t M, std::size_t N, std::size_t O>
  void matrix_test_ops_inhomogeneous_multipliable_matrices() {
    // Test NxM and MxO matrix multiplication
    {
      matrix<N, M> m1;
      matrix<M, O> m2;

      for (std::size_t i = 0; i < N; ++i) {
        for (std::size_t j = 0; j < M; ++j) {
          detray::getter::element(m1, i, j) = static_cast<scalar>(i * N + j);
        }
      }

      for (std::size_t i = 0; i < M; ++i) {
        for (std::size_t j = 0; j < O; ++j) {
          detray::getter::element(m2, i, j) = static_cast<scalar>(i * M + j);
        }
      }

      {
        matrix<N, O> r1 = m1 * m2;
        matrix<N, O> r2 =
            detray::matrix::transposed_product<false, false>(m1, m2);

        for (std::size_t i = 0; i < N; ++i) {
          for (std::size_t j = 0; j < O; ++j) {
            ASSERT_NEAR(detray::getter::element(r1, i, j),
                        detray::getter::element(r2, i, j), this->epsilon());
          }
        }
      }
    }

    // Test NxM and (OxM)^T matrix multiplication
    {
      matrix<N, M> m1;
      matrix<O, M> m2;

      for (std::size_t i = 0; i < N; ++i) {
        for (std::size_t j = 0; j < M; ++j) {
          detray::getter::element(m1, i, j) = static_cast<scalar>(i * N + j);
        }
      }

      for (std::size_t i = 0; i < O; ++i) {
        for (std::size_t j = 0; j < M; ++j) {
          detray::getter::element(m2, i, j) = static_cast<scalar>(i * O + j);
        }
      }

      {
        matrix<N, O> r1 = m1 * detray::matrix::transpose(m2);
        matrix<N, O> r2 =
            detray::matrix::transposed_product<false, true>(m1, m2);

        for (std::size_t i = 0; i < N; ++i) {
          for (std::size_t j = 0; j < O; ++j) {
            ASSERT_NEAR(detray::getter::element(r1, i, j),
                        detray::getter::element(r2, i, j), this->epsilon());
          }
        }
      }
    }

    // Test (MxN)^T and MxO matrix multiplication
    {
      matrix<M, N> m1;
      matrix<M, O> m2;

      for (std::size_t i = 0; i < M; ++i) {
        for (std::size_t j = 0; j < N; ++j) {
          detray::getter::element(m1, i, j) = static_cast<scalar>(i * M + j);
        }
      }

      for (std::size_t i = 0; i < M; ++i) {
        for (std::size_t j = 0; j < O; ++j) {
          detray::getter::element(m2, i, j) = static_cast<scalar>(i * M + j);
        }
      }

      {
        matrix<N, O> r1 = detray::matrix::transpose(m1) * m2;
        matrix<N, O> r2 =
            detray::matrix::transposed_product<true, false>(m1, m2);

        for (std::size_t i = 0; i < N; ++i) {
          for (std::size_t j = 0; j < O; ++j) {
            ASSERT_NEAR(detray::getter::element(r1, i, j),
                        detray::getter::element(r2, i, j), this->epsilon());
          }
        }
      }
    }

    // Test (MxN)^T and (OxM)^T matrix multiplication
    {
      matrix<M, N> m1;
      matrix<O, M> m2;

      for (std::size_t i = 0; i < M; ++i) {
        for (std::size_t j = 0; j < N; ++j) {
          detray::getter::element(m1, i, j) = static_cast<scalar>(i * M + j);
        }
      }

      for (std::size_t i = 0; i < O; ++i) {
        for (std::size_t j = 0; j < M; ++j) {
          detray::getter::element(m2, i, j) = static_cast<scalar>(i * O + j);
        }
      }

      {
        matrix<N, O> r1 =
            detray::matrix::transpose(m1) * detray::matrix::transpose(m2);
        matrix<N, O> r2 =
            detray::matrix::transposed_product<true, true>(m1, m2);

        for (std::size_t i = 0; i < N; ++i) {
          for (std::size_t j = 0; j < O; ++j) {
            ASSERT_NEAR(detray::getter::element(r1, i, j),
                        detray::getter::element(r2, i, j), this->epsilon());
          }
        }
      }
    }
  }
#endif
};

}  // namespace detray::test

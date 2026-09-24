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
#include "detray/definitions/detail/qualifiers.hpp"

// System include(s).
#include <cstddef>
#include <type_traits>

namespace detray::algebra::array {

/// @name Operators on 2-element arrays
/// @{

template <concepts::scalar scalar_t>
DETRAY_HOST_DEVICE constexpr std::array<scalar_t, 2> operator*(
    const std::array<scalar_t, 2> &a, scalar_t s) {
  return {a[0] * s, a[1] * s};
}

template <concepts::scalar scalar_t>
DETRAY_HOST_DEVICE constexpr std::array<scalar_t, 2> operator*(
    scalar_t s, const std::array<scalar_t, 2> &a) {
  return {s * a[0], s * a[1]};
}

template <concepts::scalar scalar_t>
DETRAY_HOST_DEVICE constexpr std::array<scalar_t, 2> operator-(
    const std::array<scalar_t, 2> &a, const std::array<scalar_t, 2> &b) {
  return {a[0] - b[0], a[1] - b[1]};
}

template <concepts::scalar scalar_t>
DETRAY_HOST_DEVICE constexpr std::array<scalar_t, 2> operator+(
    const std::array<scalar_t, 2> &a, const std::array<scalar_t, 2> &b) {
  return {a[0] + b[0], a[1] + b[1]};
}

/// @}

/// @name Operators on 3-element arrays
/// @{

template <concepts::scalar scalar_t>
DETRAY_HOST_DEVICE constexpr std::array<scalar_t, 3> operator*(
    const std::array<scalar_t, 3> &a, scalar_t s) {
  return {a[0] * s, a[1] * s, a[2] * s};
}

template <concepts::scalar scalar_t>
DETRAY_HOST_DEVICE constexpr std::array<scalar_t, 3> operator*(
    scalar_t s, const std::array<scalar_t, 3> &a) {
  return {s * a[0], s * a[1], s * a[2]};
}

template <concepts::scalar scalar_t>
DETRAY_HOST_DEVICE constexpr std::array<scalar_t, 3> operator-(
    const std::array<scalar_t, 3> &a, const std::array<scalar_t, 3> &b) {
  return {a[0] - b[0], a[1] - b[1], a[2] - b[2]};
}

template <concepts::scalar scalar_t>
DETRAY_HOST_DEVICE constexpr std::array<scalar_t, 3> operator+(
    const std::array<scalar_t, 3> &a, const std::array<scalar_t, 3> &b) {
  return {a[0] + b[0], a[1] + b[1], a[2] + b[2]};
}

/// @}

/// @name Operators on matrix
/// @{

template <typename index_t, concepts::scalar scalar_t, index_t ROWS,
          index_t COLS>
DETRAY_HOST_DEVICE constexpr std::array<std::array<scalar_t, ROWS>, COLS>
operator*(const std::array<std::array<scalar_t, ROWS>, COLS> &a, scalar_t s) {
  std::array<std::array<scalar_t, ROWS>, COLS> ret{};

  for (index_t j = 0; j < COLS; ++j) {
    for (index_t i = 0; i < ROWS; ++i) {
      ret[j][i] = a[j][i] * s;
    }
  }

  return ret;
}

template <typename index_t, concepts::scalar scalar_t, index_t ROWS,
          index_t COLS>
DETRAY_HOST_DEVICE constexpr std::array<std::array<scalar_t, ROWS>, COLS>
operator*(scalar_t s, const std::array<std::array<scalar_t, ROWS>, COLS> &a) {
  std::array<std::array<scalar_t, ROWS>, COLS> ret{};

  for (index_t j = 0; j < COLS; ++j) {
    for (index_t i = 0; i < ROWS; ++i) {
      ret[j][i] = a[j][i] * s;
    }
  }

  return ret;
}

template <typename index_t, concepts::scalar scalar_t, index_t M, index_t N,
          index_t O>
DETRAY_HOST_DEVICE constexpr std::array<std::array<scalar_t, M>, O> operator*(
    const std::array<std::array<scalar_t, M>, N> &A,
    const std::array<std::array<scalar_t, N>, O> &B) {
  std::array<std::array<scalar_t, M>, O> C{};

  for (index_t j = 0; j < O; ++j) {
    for (index_t i = 0; i < M; ++i) {
      C[j][i] = 0.f;
    }
  }

  for (index_t i = 0; i < N; ++i) {
    for (index_t j = 0; j < O; ++j) {
      for (index_t k = 0; k < M; ++k) {
        C[j][k] += A[i][k] * B[j][i];
      }
    }
  }

  return C;
}

template <typename index_t, concepts::scalar scalar_t, index_t ROWS,
          index_t COLS>
DETRAY_HOST_DEVICE constexpr std::array<std::array<scalar_t, ROWS>, COLS>
operator+(const std::array<std::array<scalar_t, ROWS>, COLS> &A,
          const std::array<std::array<scalar_t, ROWS>, COLS> &B) {
  std::array<std::array<scalar_t, ROWS>, COLS> C{};

  for (index_t j = 0; j < COLS; ++j) {
    for (index_t i = 0; i < ROWS; ++i) {
      C[j][i] = A[j][i] + B[j][i];
    }
  }

  return C;
}

template <typename index_t, concepts::scalar scalar_t, index_t ROWS,
          index_t COLS>
DETRAY_HOST_DEVICE constexpr std::array<std::array<scalar_t, ROWS>, COLS>
operator-(const std::array<std::array<scalar_t, ROWS>, COLS> &A,
          const std::array<std::array<scalar_t, ROWS>, COLS> &B) {
  std::array<std::array<scalar_t, ROWS>, COLS> C{};

  for (index_t j = 0; j < COLS; ++j) {
    for (index_t i = 0; i < ROWS; ++i) {
      C[j][i] = A[j][i] - B[j][i];
    }
  }

  return C;
}

/// @}

/// @name Operators on matrix * vector
/// @{

template <typename index_t, concepts::scalar scalar_t, index_t ROWS,
          index_t COLS>
DETRAY_HOST_DEVICE constexpr std::array<scalar_t, ROWS> operator*(
    const std::array<std::array<scalar_t, ROWS>, COLS> &a,
    const std::array<scalar_t, COLS> &b) {
  std::array<scalar_t, ROWS> ret{0};

  for (index_t j = 0; j < COLS; ++j) {
    for (index_t i = 0; i < ROWS; ++i) {
      ret[i] += a[j][i] * b[j];
    }
  }

  return ret;
}

/// @}

}  // namespace detray::algebra::array

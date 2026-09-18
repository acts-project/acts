// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// Project include(s).
#include "algebra/impl/array_soa_casts.hpp"
#include "algebra/impl/array_soa_concepts.hpp"
#include "algebra/impl/array_soa_getter.hpp"
#include "algebra/impl/array_soa_simd.hpp"
#include "algebra/impl/array_soa_transform3.hpp"
#include "detray/algebra/common/matrix.hpp"
#include "detray/algebra/common/vector.hpp"
#include "detray/algebra/concepts.hpp"
#include "detray/algebra/type_traits.hpp"

// System include(s).
#include <array>
#include <cstddef>

namespace detray {

namespace algebra::array_soa {

/// Size type for the std::array SoA storage model
using index_type = std::size_t;
/// Value type in a linear algebra vector: SoA layout
template <concepts::value T>
using value_type = T;
/// Scalar type in a linear algebra vector: SoA layout
template <concepts::value T, std::size_t W>
using scalar_type = simd<T, W>;
/// Array type used to store lane bundles or matrix columns
template <concepts::simd_scalar T, index_type N>
using storage_type = std::array<T, N>;
/// Vector type used in the std::array SoA storage model
template <concepts::value T, std::size_t W, std::size_t N>
using vector_type = algebra::storage::vector<N, simd<T, W>, storage_type>;
/// Matrix type used in the std::array SoA storage model
template <concepts::value T, std::size_t W, index_type ROWS, index_type COLS>
using matrix_type =
    algebra::storage::matrix<storage_type, simd<T, W>, ROWS, COLS>;

/// 2-element "vector" type, using a lane bundle in every element
template <concepts::value T, std::size_t W>
using vector2 = vector_type<T, W, 2>;
/// Point in 2D space, using a lane bundle in every element
template <concepts::value T, std::size_t W>
using point2 = vector2<T, W>;
/// 3-element "vector" type, using a lane bundle in every element
template <concepts::value T, std::size_t W>
using vector3 = vector_type<T, W, 3>;
/// Point in 3D space, using a lane bundle in every element
template <concepts::value T, std::size_t W>
using point3 = vector3<T, W>;
/// 6-element "vector" type, using a lane bundle in every element
template <concepts::value T, std::size_t W>
using vector6 = vector_type<T, W, 6>;
/// 8-element "vector" type, using a lane bundle in every element
template <concepts::value T, std::size_t W>
using vector8 = vector_type<T, W, 8>;

/// Geometry transformation implementation using a lane bundle in every element
template <concepts::value T, std::size_t W>
using transform3 = array_soa::math::transform3<storage_type, simd<T, W>>;

/// Element Getter
using element_getter = algebra::storage::element_getter;
/// Block Getter
using block_getter = algebra::storage::block_getter;

}  // namespace algebra::array_soa

// The default type trait macro cannot express the lane count, so the
// specialisations are written out here
namespace traits {

/// Index type
/// @{
template <concepts::value T, std::size_t W, auto N>
struct index<algebra::array_soa::vector_type<T, W, N>> {
  using type = algebra::array_soa::index_type;
};

template <concepts::value T, std::size_t W, auto ROWS, auto COLS>
struct index<algebra::array_soa::matrix_type<T, W, ROWS, COLS>> {
  using type = algebra::array_soa::index_type;
};
/// @}

/// Dimensions
/// @{
template <concepts::value T, std::size_t W, auto N>
struct dimensions<algebra::array_soa::vector_type<T, W, N>> {
  using index_type = algebra::array_soa::index_type;

  static constexpr index_type _dim{1};
  static constexpr index_type _rows{N};
  static constexpr index_type _columns{1};
};

template <concepts::value T, std::size_t W, auto ROWS, auto COLS>
struct dimensions<algebra::array_soa::matrix_type<T, W, ROWS, COLS>> {
  using index_type = algebra::array_soa::index_type;

  static constexpr index_type _dim{2};
  static constexpr index_type _rows{ROWS};
  static constexpr index_type _columns{COLS};
};
/// @}

/// Value type (single precision value) and scalar type (lane bundle)
/// @{
template <concepts::value T, std::size_t W, auto N>
struct value<algebra::array_soa::vector_type<T, W, N>> {
  using type = T;
};

template <concepts::value T, std::size_t W, auto ROWS, auto COLS>
struct value<algebra::array_soa::matrix_type<T, W, ROWS, COLS>> {
  using type = T;
};

template <concepts::value T, std::size_t W>
struct value<algebra::array_soa::simd<T, W>> {
  using type = T;
};

template <concepts::value T, std::size_t W, auto N>
struct scalar<algebra::array_soa::vector_type<T, W, N>> {
  using type = algebra::array_soa::simd<T, W>;
};

template <concepts::value T, std::size_t W, auto ROWS, auto COLS>
struct scalar<algebra::array_soa::matrix_type<T, W, ROWS, COLS>> {
  using type = algebra::array_soa::simd<T, W>;
};
/// @}

/// Compatible vector and matrix types
/// @{
template <concepts::value T, std::size_t W, auto N>
struct vector<algebra::array_soa::vector_type<T, W, N>> {
  template <typename other_T, auto other_N>
  using other_type = algebra::array_soa::vector_type<other_T, W, other_N>;

  using type = other_type<T, N>;
};

template <concepts::value T, std::size_t W, auto ROWS, auto COLS>
struct vector<algebra::array_soa::matrix_type<T, W, ROWS, COLS>> {
  template <typename other_T, auto other_N>
  using other_type = algebra::array_soa::vector_type<other_T, W, other_N>;

  using type = other_type<T, ROWS>;
};

template <concepts::value T, std::size_t W, auto ROWS, auto COLS>
struct matrix<algebra::array_soa::matrix_type<T, W, ROWS, COLS>> {
  template <typename other_T, auto other_ROWS, auto other_COLS>
  using other_type =
      algebra::array_soa::matrix_type<other_T, W, other_ROWS, other_COLS>;

  using type = algebra::array_soa::matrix_type<T, W, ROWS, COLS>;
};

template <concepts::value T, std::size_t W, auto N>
struct matrix<algebra::array_soa::vector_type<T, W, N>> {
  template <typename other_T, auto other_ROWS, auto other_COLS>
  using other_type =
      algebra::array_soa::matrix_type<other_T, W, other_ROWS, other_COLS>;

  using type = other_type<T, N, 1>;
};
/// @}

/// Getters
/// @{
template <concepts::value T, std::size_t W, auto N>
struct element_getter<algebra::array_soa::vector_type<T, W, N>> {
  using type = algebra::array_soa::element_getter;
};

template <concepts::value T, std::size_t W, auto ROWS, auto COLS>
struct element_getter<algebra::array_soa::matrix_type<T, W, ROWS, COLS>> {
  using type = algebra::array_soa::element_getter;
};

template <concepts::value T, std::size_t W, auto ROWS, auto COLS>
struct block_getter<algebra::array_soa::matrix_type<T, W, ROWS, COLS>> {
  using type = algebra::array_soa::block_getter;
};
/// @}

// Vector and storage types are different: the array of lane bundles that the
// storage operators return must also count as a vector
template <concepts::array_soa_simd T, auto N>
struct index<algebra::array_soa::storage_type<T, N>> {
  using type = algebra::array_soa::index_type;
};

template <concepts::array_soa_simd T, auto N>
struct dimensions<algebra::array_soa::storage_type<T, N>> {
  using index_type = algebra::array_soa::index_type;

  static constexpr index_type _dim{1};
  static constexpr index_type _rows{N};
  static constexpr index_type _columns{1};
};

}  // namespace traits

}  // namespace detray

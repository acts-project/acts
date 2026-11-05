// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// Project include(s).
#include "algebra/impl/vc_aos_approximately_equal.hpp"
#include "algebra/impl/vc_aos_concepts.hpp"
#include "algebra/impl/vc_aos_getter.hpp"
#include "algebra/impl/vc_aos_transform3.hpp"
#include "detray/algebra/common/matrix.hpp"
#include "detray/algebra/common/vector.hpp"
#include "detray/algebra/concepts.hpp"
#include "detray/algebra/type_traits.hpp"

// System include(s).
#include <array>
#include <cstddef>

// Vc include(s).
#ifdef _MSC_VER
#pragma warning(push, 0)
#endif  // MSVC
#include <Vc/Vc>
#ifdef _MSC_VER
#pragma warning(pop)
#endif  // MSVC

namespace detray {

namespace algebra::vc_aos {

/// Size type for Vc storage model
using index_type = std::size_t;
/// Value type for Vc storage model
template <concepts::value T>
using value_type = T;
/// Scalar type for Vc storage model
template <concepts::value T>
using scalar_type = T;
/// Array type used to store Vc::Vectors
template <concepts::value T, index_type N>
using storage_type = Vc::SimdArray<T, N>;
/// Value type in a linear algebra vector: AoS layout
template <concepts::value T>
using value_type = T;
/// Scalar type in a linear algebra vector: AoS layout
template <concepts::value T>
using scalar_type = T;
/// Vector type used in the Vc AoS storage model
template <concepts::value T, std::size_t N>
using vector_type = algebra::storage::vector<N, T, storage_type>;
/// Matrix type used in the Vc AoS storage model
template <concepts::value T, index_type ROWS, index_type COLS>
using matrix_type = algebra::storage::matrix<storage_type, T, ROWS, COLS>;

/// 2-element "vector" type, using @c Vc::SimdArray
template <concepts::value T>
using vector2 = vector_type<T, 2>;
/// Point in 2D space, using @c Vc::SimdArray
template <concepts::value T>
using point2 = vector2<T>;
/// 3-element "vector" type, using @c Vc::SimdArray
template <concepts::value T>
using vector3 = vector_type<T, 3>;
/// Point in 3D space, using @c Vc::SimdArray
template <concepts::value T>
using point3 = vector3<T>;
/// 6-element "vector" type, using @c Vc::SimdArray
template <concepts::value T>
using vector6 = vector_type<T, 6>;
/// 8-element "vector" type, using @c Vc::SimdArray
template <concepts::value T>
using vector8 = vector_type<T, 8>;

/// Geometry transformation implementation using @c Vc::SimdArray
template <concepts::value T>
using transform3 = vc_aos::math::transform3<vc_aos::storage_type, T>;

/// Element Getter
using element_getter = algebra::storage::element_getter;
/// Block Getter
using block_getter = algebra::storage::block_getter;

}  // namespace algebra::vc_aos

DETRAY_ALGEBRA_DEFINE_TYPE_TRAITS(algebra::vc_aos)

namespace traits {

/// Index type associated with a Vc AoS storage type.
template <typename T, auto N>
struct index<algebra::vc_aos::storage_type<T, N>> {
  /// Index type used by the Vc AoS backend.
  using type = algebra::vc_aos::index_type;
};

/// Index trait for a fixed-size Vc vector.
template <typename T, auto N>
struct index<Vc_1::Vector<T, Vc_1::simd_abi::fixed_size<N>>> {
  /// Index type used by the Vc AoS backend.
  using type = algebra::vc_aos::index_type;
};

/// Scalar type associated with a Vc AoS storage type.
template <typename T, auto N>
struct value<algebra::vc_aos::storage_type<T, N>> {
  /// Scalar element type carried by the Vc array.
  using type = T;
};

/// Value trait for a fixed-size Vc vector.
template <typename T, auto N>
struct value<Vc_1::Vector<T, Vc_1::simd_abi::fixed_size<N>>> {
  /// Scalar element type carried by the Vc vector.
  using type = T;
};

// Vector and storage types are different
/// Dimensions of Vc AoS storage, represented as one column of `N` rows.
template <typename T, auto N>
struct dimensions<algebra::vc_aos::storage_type<T, N>> {
  /// Index type used for the storage dimensions.
  using index_type = algebra::vc_aos::index_type;

  /// Number of logical dimensions represented by the storage.
  static constexpr index_type _dim{1};
  /// Number of stored vector components.
  static constexpr index_type _rows{N};
  /// Number of storage columns.
  static constexpr index_type _columns{1};
};

// Vector and storage types are different
/// Dimensions of a fixed-size Vc vector, represented as one column.
template <typename T, auto N>
struct dimensions<Vc_1::Vector<T, Vc_1::simd_abi::fixed_size<N>>> {
  /// Index type used for the vector dimensions.
  using index_type = algebra::vc_aos::index_type;

  /// Number of logical dimensions represented by the vector.
  static constexpr index_type _dim{1};
  /// Number of vector components.
  static constexpr index_type _rows{N};
  /// Number of vector columns.
  static constexpr index_type _columns{1};
};

}  // namespace traits

}  // namespace detray

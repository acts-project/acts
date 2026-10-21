// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// Project include(s)
#include "algebra/impl/array_getter.hpp"
#include "detray/algebra/concepts.hpp"
#include "detray/algebra/generic/generic.hpp"
#include "detray/algebra/type_traits.hpp"

// System include(s).
#include <array>
#include <cstddef>

namespace detray {

namespace algebra::array {

/// size type for Array storage model
using index_type = std::size_t;
/// Value type for Array storage model
template <concepts::value T>
using value_type = T;
/// Scalar type for Array storage model
template <concepts::value T>
using scalar_type = T;
/// Array type used in the Array storage model
template <typename T, index_type N>
using storage_type = std::array<T, N>;
/// Vector type used in the Array storage model
template <concepts::scalar T, std::size_t N>
using vector_type = storage_type<T, N>;
/// Matrix type used in the Array storage model
template <concepts::scalar T, index_type ROWS, index_type COLS>
using matrix_type = storage_type<storage_type<T, ROWS>, COLS>;

/// 3-element "vector" type, using @c std::array
template <concepts::scalar T>
using vector3 = storage_type<T, 3>;
/// Point in 3D space, using @c std::array
template <concepts::scalar T>
using point3 = vector3<T>;
/// 2-element "vector" type, using @c std::array
template <concepts::scalar T>
using vector2 = storage_type<T, 2>;
/// Point in 2D space, using @c std::array
template <concepts::scalar T>
using point2 = vector2<T>;

/// Generic geometry transformation implementation using @c std::array
template <concepts::scalar T>
using transform3 =
    algebra::generic::math::transform3<algebra::array::index_type, T,
                                       algebra::array::matrix_type,
                                       algebra::array::storage_type>;

/// Element Getter
using element_getter = array::storage::element_getter;
/// Block Getter
using block_getter = array::storage::block_getter;

}  // namespace algebra::array

DETRAY_ALGEBRA_DEFINE_TYPE_TRAITS(algebra::array)

}  // namespace detray

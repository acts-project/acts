// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "algebra/impl/detail/fastor_matrix_wrapper.hpp"
#include "algebra/impl/detail/fastor_vector_wrapper.hpp"
#include "algebra/impl/fastor_getter.hpp"
#include "algebra/impl/fastor_transform3.hpp"
#include "detray/algebra/concepts.hpp"
#include "detray/algebra/type_traits.hpp"

// System include(s).
#include <cstddef>

namespace detray {

namespace algebra::fastor {

/// size type for Fastor storage model
using index_type = std::size_t;
/// Value type for Fastor storage model
template <concepts::value T>
using value_type = T;
/// Scalar type for Fastor storage model
template <concepts::value T>
using scalar_type = T;
/// Array type used in the Fastor storage model
template <concepts::scalar T, index_type N>
using storage_type = fastor::Vector<T, N>;
/// Vector type used in the Fastor storage model
template <concepts::scalar T, index_type N>
using vector_type = storage_type<T, N>;
/// Matrix type used in the Fastor storage model
template <concepts::scalar T, index_type ROWS, index_type COLS>
using matrix_type = fastor::Matrix<T, ROWS, COLS>;

/// 3-element "vector" type, using @c Fastor::Tensor
template <concepts::scalar T>
using vector3 = storage_type<T, 3>;
/// Point in 3D space, using @c Fastor::Tensor
template <concepts::scalar T>
using point3 = vector3<T>;
/// 2-element "vector" type, using @c Fastor::Tensor
template <concepts::scalar T>
using vector2 = storage_type<T, 2>;
/// Point in 2D space, using @c Fastor::Tensor
template <concepts::scalar T>
using point2 = vector2<T>;

/// Geometry transformation implementation using @c Fastor::Tensor
template <concepts::scalar T>
using transform3 = fastor::math::transform3<T>;

/// Element Getter
using element_getter = fastor::storage::element_getter;
/// Block Getter
using block_getter = fastor::storage::block_getter;

}  // namespace algebra::fastor

DETRAY_ALGEBRA_DEFINE_TYPE_TRAITS(algebra::fastor)

}  // namespace detray

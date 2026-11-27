// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// Project include(s).
#include "detray/algebra/common/storage/vector.hpp"
#include "detray/algebra/concepts.hpp"
#include "detray/algebra/type_traits.hpp"

// Vc include(s).
#ifdef _MSC_VER
#pragma warning(push, 0)
#endif  // MSVC
#include <Vc/Vc>
#ifdef _MSC_VER
#pragma warning(pop)
#endif  // MSVC

// System include(s).
#include <concepts>

namespace detray::concepts {

/// Vc SIMD types
template <typename T>
concept vc_simd_vector = Vc::is_simd_vector<T>::value;

/// Vc AoS storage vector: a storage vector with a non-SIMD scalar
///
/// @note The SoA plugins also use storage vectors, but with a SIMD scalar.
/// These do not match, so that the Vc AoS overloads do not clash with the
/// SoA overloads in the same translation unit.
template <typename T>
concept simd_storage_vector = algebra::detail::is_storage_vector_v<T> &&
                              !concepts::simd_scalar<traits::scalar_t<T>>;

/// Vc AoS vector
template <typename T>
concept vc_aos_vector = (vc_simd_vector<T> || simd_storage_vector<T>);

}  // namespace detray::concepts

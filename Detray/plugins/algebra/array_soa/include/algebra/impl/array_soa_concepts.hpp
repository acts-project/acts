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

// System include(s).
#include <cstddef>
#include <type_traits>

namespace detray {

namespace algebra::array_soa {

/// Lane bundle: one value per SIMD lane, stored in a @c std::array
template <concepts::value T, std::size_t W>
struct simd;

/// Result of a lane-wise comparison
template <std::size_t W>
struct mask;

namespace detail {

template <typename T>
struct is_simd : public std::false_type {};

template <concepts::value T, std::size_t W>
struct is_simd<simd<T, W>> : public std::true_type {};

template <typename T>
struct is_mask : public std::false_type {};

template <std::size_t W>
struct is_mask<mask<W>> : public std::true_type {};

}  // namespace detail

}  // namespace algebra::array_soa

namespace concepts {

/// std::array based SIMD lane bundle
template <typename T>
concept array_soa_simd =
    algebra::array_soa::detail::is_simd<std::remove_cvref_t<T>>::value;

/// std::array based SIMD mask
template <typename T>
concept array_soa_mask =
    algebra::array_soa::detail::is_mask<std::remove_cvref_t<T>>::value;

}  // namespace concepts

}  // namespace detray

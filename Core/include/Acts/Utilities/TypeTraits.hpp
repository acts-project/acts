// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include <type_traits>
#include <variant>

namespace Acts {
/// @brief Type trait that adds const qualifier to a type based on a boolean condition
/// @tparam C Boolean condition determining whether to add const
/// @tparam T Type to potentially make const
template <bool C, typename T>
using const_if_t = std::conditional_t<C, const T, T>;

/// @brief Type trait that determines whether a type can be used as an
/// alternative in a std::variant.
/// @tparam T The type to query
/// @tparam Variant The variant type to inspect
template <typename T, typename Variant>
struct is_variant_alternative : std::false_type {};

template <typename T, typename... Types>
struct is_variant_alternative<T, std::variant<Types...>>
    : std::bool_constant<(
          std::is_same_v<std::remove_cv_t<std::remove_reference_t<T>>,
                         std::remove_cv_t<std::remove_reference_t<Types>>>
          || ...)> {};

/// @brief Concept that is satisfied if a type can be used as an alternative in
/// a std::variant.
template <typename T, typename Variant>
concept VariantCompatible =
    is_variant_alternative<T, std::remove_cv_t<std::remove_reference_t<Variant>>>::value;
}  // namespace Acts

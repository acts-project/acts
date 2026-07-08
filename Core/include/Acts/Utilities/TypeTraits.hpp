// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Utilities/PointerTraits.hpp"

#include <type_traits>
#include <variant>
namespace Acts {
/// Type trait that adds const qualifier to a type based on a boolean condition
/// @tparam C Boolean condition determining whether to add const
/// @tparam T Type to potentially make const
template <bool C, typename T>
using const_if_t = std::conditional_t<C, const T, T>;

/// Helper type traits to specify whether an object of type Obj_t can be packed
/// into a std::variant<T, T1, ...>. There we, distinguish two cases:
///         1) The types T, T1 are ordinary data types -> T or T1 must be of the
///            same type as Obj_t
///         2) The types T, T1 are pointer types --> Allow that the data type
///            to which Obj_t points inherits from the data type pointed
///            to by T or by T1
namespace detail {

/// Define the castable type trait between two template parameters
//  to be the same if their const pendant is the same.
/// @tparam T First data type to compare
/// @tparam T1 Second data type to compare
template <typename T, typename T1>
struct is_castable : std::is_same<std::add_const_t<T>, std::add_const_t<T1>> {};

/// Specification if the types T and T1 actually are pointers. The underlying
/// pointer data type of T must be a base class of T1
/// @tparam T First pointer data type to compare
/// @tparam T1 Second pointer data type to compare
template <PointerConcept T, PointerConcept T1>
struct is_castable<T, T1>
    : std::is_base_of<std::add_const_t<RemovePointer_t<T>>,
                      std::add_const_t<RemovePointer_t<T1>>> {};
/// Check that T is castable to at least one of the passed type parameters
/// @tparam T The template parameter to check
/// @tparam Args The list of available template parameters out of which one is
///              supposed to satisfy the is_castable comparison with T
template <typename T, typename... Args>
struct holds_args : std::bool_constant<(is_castable<Args, T>::value || ...)> {};

/// Baseline type trait to check whether the template parameter T is amongst the
/// variant arguments. General case to indicate false
/// @tparam T The template argument to check
/// @tparam Variants Any arbitrary data type
template <typename T, typename Variants>
struct is_variant_alternative : std::false_type {};
/// Type trait specification for std::variant. Check whether T is compatible
/// with one of the arguments in the variant or whether it is the variant type
/// itself.
/// @tparam T The template argument to check
/// @tparam Variants The list of available variant arguments
template <typename T, typename... Variants>
struct is_variant_alternative<T, std::variant<Variants...>>
    : std::bool_constant<holds_args<T, Variants...>::value ||
                         std::is_same_v<T, std::variant<Variants...>>> {};

}  // namespace detail

/// Concept whether the template parameter T is part of a
/// std variant
/// @tparam Variant std::variant template parameter
/// @tparam  T The template parameter to check
template <typename Variant, typename T>
concept isVariantCompatible = detail::is_variant_alternative<T, Variant>::value;

}  // namespace Acts

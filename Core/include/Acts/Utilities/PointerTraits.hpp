// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include <concepts>
#include <type_traits>

namespace Acts {

/// @brief The Pointer concept is an extension of the usual std::is_pointer_v type trait to
///         also include the smart pointers like
///     std::shared_ptr<T>,
///     std::unique_ptr<T>
///  The smart pointer is required to have an element_type typedef indicating
/// over which data type the pointer is constructed, the arrow operator
///
///      T* operator->() const;
///  and also the dereference operator
///
///      T& operator*() const
template <typename Pointer_t>
concept SmartPointerConcept = requires(Pointer_t ptr) {
  typename Pointer_t::element_type;
  /// @brief arrow operator element_type* operator->() const;
  {
    ptr.operator->()
  } -> std::same_as<std::add_pointer_t<typename Pointer_t::element_type>>;
  /// @brief dereference operator element_type& operator*() const;
  { ptr.operator*() } -> std::same_as<typename Pointer_t::element_type&>;
  /// @brief Simple cast to check for if(ptr)
  { ptr.operator bool() };
};
template <typename Pointer_t>
concept PointerConcept =
    (std::is_pointer_v<Pointer_t> || SmartPointerConcept<Pointer_t>);
/// @brief Introduce the Acts version of the pointer remove type trait because we want to
///        fetch the underlying type for the pointer concept and std::library
///        does not allow for an extension of the std::remove_pointer;
template <typename T>
struct RemovePointer {
  /// Type alias for the original type (non-pointer case)
  using type = T;
};

/// @brief  This specialization allows std::remove_pointer to work with types satisfying
///         Acts::SmartPointerConcept, similar to how it works with raw pointers
template <SmartPointerConcept T>
struct RemovePointer<T> {
  /// Type alias for the element type pointed to by smart pointer
  using type = typename T::element_type;
};
/// @brief ordinary specialization for pointers
template <PointerConcept T>
struct RemovePointer<T> {
  /// Type alias for the type pointed to by raw pointer
  using type = std::remove_pointer_t<T>;
};
/// Helper type alias for removing pointer from type
template <typename T>
using RemovePointer_t = RemovePointer<T>::type;

/// @brief Construct a const pointer dereference
template <PointerConcept T>
using ConstDeRef_t = std::add_const_t<RemovePointer_t<T>>&;

}  // namespace Acts

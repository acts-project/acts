// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Utilities/PointerTraits.hpp"

#include <concepts>
#include <type_traits>
#include <utility>

namespace Acts {

/// @brief Non-owning holder for referencing a backend.
/// @tparam T Backend type.
///
/// The referenced backend must outlive the holder.
template <typename T>
struct RefHolder {
  using element_type = T;

  T* ptr;

  explicit RefHolder(T* _ptr) : ptr{_ptr} {}
  explicit RefHolder(T& ref) : ptr{&ref} {}

  const T& operator*() const { return *ptr; }
  T& operator*() { return *ptr; }

  const T* operator->() const { return ptr; }
  T* operator->() { return ptr; }

  explicit operator bool() const { return ptr != nullptr; }
};

/// @brief Non-owning holder for referencing a backend with const access.
/// @tparam T Backend type.
///
/// The referenced backend must outlive the holder.
template <typename T>
struct ConstRefHolder {
  using element_type = std::add_const_t<T>;

  const T* ptr;

  explicit ConstRefHolder(const T* _ptr) : ptr{_ptr} {}
  explicit ConstRefHolder(const T& ref) : ptr{&ref} {}

  const T& operator*() const { return *ptr; }

  const T* operator->() const { return ptr; }

  explicit operator bool() const { return ptr != nullptr; }
};

/// @brief Owning holder that stores a backend by value.
/// @tparam T Backend type.
///
/// The backend is moved into the holder and owned for its lifetime.
template <typename T>
struct ValueHolder {
  using element_type = T;

  T val;

  // Let's be clear with the user that we take the ownership
  // Only require rvalues and avoid hidden copies
  ValueHolder(T& _val) = delete;
  // @FIXME: Ideally we want this to be explicit, but cannot be explicit,
  // because using an explicit constructor and a deduction guide leads to
  // a SEGFAULT in GCC11 (an up?). Re-evaluate down the line
  /* explicit */ ValueHolder(T&& _val) : val{std::move(_val)} {}  // NOLINT

  // Does it make sense to allow copy operations?

  const T& operator*() const { return val; }
  T& operator*() { return val; }

  const T* operator->() const { return &val; }
  T* operator->() { return &val; }

  explicit operator bool() const { return true; }
};

/// @brief Concept for holder templates that provide pointer-like access.
/// @tparam Holder Holder template to instantiate.
/// @tparam T Backend type.
///
/// This concept is satisfied by RefHolder, ConstRefHolder, ValueHolder, and
/// smart pointers such as std::shared_ptr and std::unique_ptr.
template <template <typename> class Holder, typename T>
concept HolderFor =
    std::move_constructible<Holder<T>> && PointerConcept<Holder<T>>;

}  // namespace Acts

namespace Acts::detail {
template <typename T>
using RefHolder = Acts::RefHolder<T>;

template <typename T>
using ConstRefHolder = Acts::ConstRefHolder<T>;

template <typename T>
using ValueHolder = Acts::ValueHolder<T>;
}  // namespace Acts::detail

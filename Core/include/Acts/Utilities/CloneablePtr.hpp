// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include <functional>
#include <memory>
#include <type_traits>
#include <utility>

namespace Acts {

/// @brief A copyable smart pointer that uses a cloner function to copy the
/// managed object.
///
/// This enables polymorphic value semantics: you can copy a pointer to a base
/// class by invoking a stored cloner function (typically calling a virtual
/// `clone()` method or a copy constructor).
///
/// @tparam T The type of the managed object
template <typename T>
class CloneablePtr {
 public:
  /// @brief The type of the cloner function
  using Cloner = std::function<std::unique_ptr<T>(const T&)>;

  /// @brief Default constructor, creates a null pointer
  CloneablePtr() = default;

  /// @brief Construct from a unique_ptr with a custom cloner
  /// @param ptr The unique_ptr to take ownership of
  /// @param cloner The cloner function
  CloneablePtr(std::unique_ptr<T> ptr, Cloner cloner)
      : m_ptr(std::move(ptr)), m_cloner(std::move(cloner)) {}

  /// @brief Construct from a unique_ptr using copy construction as the cloner
  /// @param ptr The unique_ptr to take ownership of
  /// @note Only available when T is copy-constructible
  explicit CloneablePtr(std::unique_ptr<T> ptr)
    requires std::is_copy_constructible_v<T>
      : m_ptr(std::move(ptr)),
        m_cloner([](const T& src) { return std::make_unique<T>(src); }) {}

  /// @brief Construct by taking ownership of a raw pointer with a custom cloner
  /// @param raw The raw pointer to take ownership of
  /// @param cloner The cloner function
  CloneablePtr(T* raw, Cloner cloner)
      : m_ptr(raw), m_cloner(std::move(cloner)) {}

  /// @brief Construct by taking ownership of a raw pointer using copy
  /// construction as the cloner
  /// @param raw The raw pointer to take ownership of
  /// @note Only available when T is copy-constructible
  explicit CloneablePtr(T* raw)
    requires std::is_copy_constructible_v<T>
      : m_ptr(raw),
        m_cloner([](const T& src) { return std::make_unique<T>(src); }) {}

  /// @brief Copy constructor. Invokes the cloner if the source is non-null.
  /// @param other The CloneablePtr to copy from
  CloneablePtr(const CloneablePtr& other)
      : m_ptr(other.m_ptr ? other.m_cloner(*other.m_ptr) : nullptr),
        m_cloner(other.m_cloner) {}

  /// @brief Copy assignment. Invokes the cloner if the source is non-null.
  /// @param other The CloneablePtr to copy from
  /// @return Reference to this
  CloneablePtr& operator=(const CloneablePtr& other) {
    if (this != &other) {
      m_ptr = other.m_ptr ? other.m_cloner(*other.m_ptr) : nullptr;
      m_cloner = other.m_cloner;
    }
    return *this;
  }

  /// @brief Move constructor
  CloneablePtr(CloneablePtr&&) = default;

  /// @brief Move assignment
  /// @return Reference to this
  CloneablePtr& operator=(CloneablePtr&&) = default;

  /// @brief Destructor
  ~CloneablePtr() = default;

  /// @brief Dereference operator
  /// @return Reference to the managed object
  T& operator*() const { return *m_ptr; }

  /// @brief Arrow operator
  /// @return Pointer to the managed object
  T* operator->() const { return m_ptr.get(); }

  /// @brief Boolean conversion, true if non-null
  explicit operator bool() const { return m_ptr != nullptr; }

  /// @brief Comparison with nullptr
  friend bool operator==(const CloneablePtr& lhs, std::nullptr_t) {
    return lhs.m_ptr == nullptr;
  }

  /// @brief Get the raw pointer
  /// @return Pointer to the managed object, or nullptr
  T* get() const { return m_ptr.get(); }

  /// @brief Release ownership of the managed object
  /// @return Pointer to the formerly managed object
  T* release() {
    m_cloner = nullptr;
    return m_ptr.release();
  }

  /// @brief Reset the managed object
  /// @param ptr The new raw pointer to manage (default nullptr)
  void reset(T* ptr = nullptr) {
    m_ptr.reset(ptr);
    if (m_ptr == nullptr) {
      m_cloner = nullptr;
    }
  }

 private:
  std::unique_ptr<T> m_ptr;
  Cloner m_cloner;
};

}  // namespace Acts

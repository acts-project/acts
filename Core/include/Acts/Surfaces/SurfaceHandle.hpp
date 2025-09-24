// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include <concepts>
#include <functional>
#include <memory>

namespace Acts {

class Surface;

/// @class MaybeSharedPtr
///
/// @brief Handle for managing Surface ownership and access
///
/// This class wraps the shared pointer semantics for Surface objects
/// while keeping the implementation details private. It provides
/// the same interface as shared_ptr but allows future flexibility
/// in changing the underlying ownership model.
///
/// @tparam T Surface type (Surface or derived classes)
template <class T>
class MaybeSharedPtr {
  /// @brief Allow conversion between different SurfaceHandle types
  template <class U>
  friend class MaybeSharedPtr;

  /// @brief Friend functions for casting
  template <class U, class V>
  friend MaybeSharedPtr<U> static_handle_cast(const MaybeSharedPtr<V>& handle);
  template <class U, class V>
  friend MaybeSharedPtr<U> dynamic_handle_cast(const MaybeSharedPtr<V>& handle);

 public:
  /// @brief Type alias for the element type
  using element_type = T;

  /// @brief Default constructor - creates empty handle
  MaybeSharedPtr() = default;

  /// @brief Construct from shared_ptr
  /// @param ptr Shared pointer to wrap
  explicit MaybeSharedPtr(std::shared_ptr<T> ptr) : m_ptr(std::move(ptr)) {}

  /// NOLINTBEGIN(google-explicit-constructor)
  MaybeSharedPtr(std::nullptr_t /*null*/) : m_ptr(nullptr) {}

  /// @brief Construct from compatible MaybeSharedPtr (e.g., derived to base, non-const to const)
  /// @param other Handle to convert from
  template <class U>
    requires std::convertible_to<U*, T*>
  MaybeSharedPtr(const MaybeSharedPtr<U>& other) : m_ptr(other.m_ptr) {}

  /// @brief Construct from compatible MaybeSharedPtr (move)
  /// @param other Handle to convert from
  template <class U>
    requires std::convertible_to<U*, T*>
  MaybeSharedPtr(MaybeSharedPtr<U>&& other) : m_ptr(std::move(other.m_ptr)) {}
  /// NOLINTEND(google-explicit-constructor)

  /// @brief Copy constructor
  MaybeSharedPtr(const MaybeSharedPtr&) = default;

  /// @brief Move constructor
  MaybeSharedPtr(MaybeSharedPtr&&) = default;

  /// @brief Copy assignment
  MaybeSharedPtr& operator=(const MaybeSharedPtr&) = default;

  /// @brief Move assignment
  MaybeSharedPtr& operator=(MaybeSharedPtr&&) = default;

  /// @brief Destructor
  ~MaybeSharedPtr() = default;

  /// @brief Dereference operator
  /// @return Reference to the managed object
  T& operator*() const { return *m_ptr; }

  /// @brief Arrow operator
  /// @return Pointer to the managed object
  T* operator->() const { return m_ptr.get(); }

  /// @brief Get raw pointer
  /// @return Raw pointer to the managed object
  T* get() const { return m_ptr.get(); }

  /// @brief Check if handle is not empty
  /// @return true if handle manages an object
  explicit operator bool() const noexcept { return static_cast<bool>(m_ptr); }

  /// @brief Reset the handle to empty state
  void reset() { m_ptr.reset(); }

  /// @brief Reset with new pointer
  /// @param ptr New pointer to manage
  void reset(std::shared_ptr<T> ptr) { m_ptr = std::move(ptr); }

  /// @brief Get weak_ptr to the managed object
  /// @return std::weak_ptr to the managed object
  std::weak_ptr<T> weak_ptr() const { return m_ptr; }

  /// @brief Swap with another handle
  /// @param other Handle to swap with
  void swap(MaybeSharedPtr& other) noexcept { m_ptr.swap(other.m_ptr); }

  /// @brief Get use count
  /// @return Number of shared owners
  long use_count() const noexcept { return m_ptr.use_count(); }

  /// @brief Check if this is the only owner
  /// @return true if use_count() == 1
  bool unique() const noexcept { return m_ptr.unique(); }

  /// @brief Equality comparison
  /// @param other Handle to compare with
  /// @return true if both handles point to the same object
  bool operator==(const MaybeSharedPtr& other) const noexcept {
    return m_ptr == other.m_ptr;
  }

  /// @brief Inequality comparison
  /// @param other Handle to compare with
  /// @return true if handles point to different objects
  bool operator!=(const MaybeSharedPtr& other) const noexcept {
    return m_ptr != other.m_ptr;
  }

  /// @brief Less-than comparison for container usage
  /// @param other Handle to compare with
  /// @return true if this handle's pointer is less than other's
  bool operator<(const MaybeSharedPtr& other) const noexcept {
    return m_ptr < other.m_ptr;
  }

  /// @brief Compare with nullptr
  /// @param ptr nullptr
  /// @return true if handle is empty
  bool operator==(std::nullptr_t) const noexcept { return !m_ptr; }

  /// @brief Compare with nullptr
  /// @param ptr nullptr
  /// @return true if handle is not empty
  bool operator!=(std::nullptr_t) const noexcept {
    return static_cast<bool>(m_ptr);
  }

 private:
  /// Internal shared pointer
  std::shared_ptr<T> m_ptr;
};

/// @brief Compare nullptr with handle
template <class T>
bool operator==(std::nullptr_t, const MaybeSharedPtr<T>& handle) noexcept {
  return handle == nullptr;
}

/// @brief Compare nullptr with handle
template <class T>
bool operator!=(std::nullptr_t, const MaybeSharedPtr<T>& handle) noexcept {
  return handle != nullptr;
}

/// @brief Swap two handles
template <class T>
void swap(MaybeSharedPtr<T>& lhs, MaybeSharedPtr<T>& rhs) noexcept {
  lhs.swap(rhs);
}

/// @brief Static cast between MaybeSharedPtr types
/// @tparam U Target type
/// @tparam T Source type
/// @param handle Source handle
/// @return MaybeSharedPtr of target type
template <class U, class T>
MaybeSharedPtr<U> static_handle_cast(const MaybeSharedPtr<T>& handle) {
  return MaybeSharedPtr<U>(std::static_pointer_cast<U>(handle.m_ptr));
}

/// @brief Dynamic cast between MaybeSharedPtr types (returns empty handle if cast fails)
/// @tparam U Target type
/// @tparam T Source type
/// @param handle Source handle
/// @return MaybeSharedPtr of target type or empty handle
template <class U, class T>
MaybeSharedPtr<U> dynamic_handle_cast(const MaybeSharedPtr<T>& handle) {
  return MaybeSharedPtr<U>(std::dynamic_pointer_cast<U>(handle.m_ptr));
}

/// Type aliases for common use cases
using SurfacePtr = MaybeSharedPtr<Surface>;
using ConstSurfacePtr = MaybeSharedPtr<const Surface>;

/// @brief Stream operator for MaybeSharedPtr
template <class T>
std::ostream& operator<<(std::ostream& os, const MaybeSharedPtr<T>& handle) {
  return os << handle.get();
}

}  // namespace Acts

/// Hash specialization for MaybeSharedPtr
namespace std {

template <class T>
struct hash<Acts::MaybeSharedPtr<T>> {
  std::size_t operator()(const Acts::MaybeSharedPtr<T>& handle) const noexcept {
    return std::hash<T*>{}(handle.get());
  }
};

}  // namespace std

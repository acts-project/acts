// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// System include(s)
#include <any>
#include <cstddef>
#include <memory>
#include <ostream>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <unordered_map>
#include <utility>

namespace detray::test {

/// A container to store arbitrary objects with ownership transfer.
///
/// This is an append-only container that takes ownership of the objects
/// added to it. Once an object has been added, it can only be read but not
/// be modified. Trying to replace an existing object is considered an error.
/// Its lifetime is bound to the lifetime of the white board.
/// @see
/// https://github.com/acts-project/acts/blob/main/Examples/Framework/include/detray::test/Framework/WhiteBoard.hpp
class whiteboard {
 public:
  whiteboard() = default;

  // A whiteboard holds unique elements and can not be copied
  whiteboard(const whiteboard& other) = delete;
  whiteboard& operator=(const whiteboard&) = delete;

  bool exists(const std::string& name) const;

  /// Store an object on the white board and transfer ownership.
  ///
  /// @param name Non-empty identifier to store it under
  /// @param object Movable reference to the transferable object
  /// @throws std::invalid_argument on empty or duplicate name
  template <typename T>
  void add(const std::string& name, T&& object);

  /// Get access to a stored object.
  ///
  /// @param[in] name Identifier for the object
  /// @return reference to the stored object
  /// @throws std::out_of_range if no object is stored under the requested
  /// name
  template <typename T>
  const T& get(const std::string& name) const;

  /// Get access to a stored object - non-const
  template <typename T>
  T& get(const std::string& name);

 private:
  /// Backend storage
  std::unordered_map<std::string, std::any> m_store{};
};

template <typename T>
inline void detray::test::whiteboard::add(const std::string& name, T&& object) {
  if (name.empty()) {
    throw std::invalid_argument("Object cannot have an empty name");
  }
  if (m_store.contains(name)) {
    throw std::invalid_argument("Object '" + name + "' already exists");
  }
  m_store.emplace(name, std::forward<T>(object));
}

template <typename T>
inline const T& detray::test::whiteboard::get(const std::string& name) const
    noexcept(false) {
  auto it = m_store.find(name);
  if (it == m_store.end()) {
    throw std::out_of_range("Object '" + name + "' does not exists");
  }
  // Try to retrieve the value as the requested type
  return std::any_cast<const T&>(it->second);
}

template <typename T>
inline T& detray::test::whiteboard::get(const std::string& name) noexcept(
    false) {
  auto it = m_store.find(name);
  if (it == m_store.end()) {
    throw std::out_of_range("Object '" + name + "' does not exists");
  }
  // Try to retrieve the value as the requested type
  return std::any_cast<T&>(it->second);
}

inline bool detray::test::whiteboard::exists(const std::string& name) const {
  return m_store.contains(name);
}

}  // namespace detray::test

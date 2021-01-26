// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <any>

namespace Acts {

/// Strong type wrapper around std::any. This has all the flexibility of
/// std::any, but it does not auto construct from anything. You have to call the
/// make static factory manually with the target type, to populate the internal
/// std::any. You can then access and modify the any as desired.
///
/// @note This is used for the context types, and should probably not used
/// outside of this use-case.
class ContextType {
 public:
  /// Default constructor, does nothing
  ///
  ContextType() {}

  /// Static factory method from arguments.
  ///
  /// @tparam T The underlying type to construct
  /// @tparam Args Types of construction arguments
  /// @param args Values of construction arguments
  /// @return ContextType The constructed context instance
  //   template <typename T, typename... Args>
  //   static ContextType make(Args&&... args) {
  //     ContextType ctx;
  //     ctx.m_data.emplace<T>(std::forward<Args>(args)...);
  //     return ctx;
  //   }

  //   template <typename T>
  //   static ContextType make(T&& value) {
  //     ContextType ctx;
  //     ctx.m_data = value;
  //     return ctx;
  //   }

  //   template <typename T>
  //   static ContextType make(const T& value) {
  //     ContextType ctx;
  //     ctx.m_data = value;
  //     return ctx;
  //   }

  template <typename T>
  explicit ContextType(T&& value) : m_data{std::move(value)} {}

  template <typename T>
  ContextType& operator=(T&& value) {
    m_data = std::move(value);
    return *this;
  }

  template <typename T>
  ContextType& operator=(const T& value) {
    m_data = std::move(value);
    return *this;
  }

  /// Mutable access to the contained any storage
  ///
  /// @return std::any&
  std::any& any() { return m_data; }

  /// Immutable access to the ocntained any storage
  ///
  /// @return const std::any&
  const std::any& any() const { return m_data; }

 private:
  std::any m_data;
};

}  // namespace Acts
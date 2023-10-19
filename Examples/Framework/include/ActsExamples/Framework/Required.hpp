// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <iostream>
#include <optional>
#include <stdexcept>

namespace ActsExamples {

/// A simple std::optional like container which can only be overwritten by a
/// defined value.
template <typename T>
class Required {
 public:
  using value_type = T;

  constexpr Required() = default;

  constexpr Required(const T& data) : m_data(data) {}
  constexpr Required(T&& data) : m_data(std::move(data)) {}

  constexpr Required(const Required& other) {
    std::cout << "copy construct" << std::endl;
    check(other.m_data);
    m_data = other.m_data;
  }
  constexpr Required(Required&& other) {
    std::cout << "move construct" << std::endl;
    check(other.m_data);
    m_data = std::move(other.m_data);
  }

  Required& operator=(const Required& other) {
    std::cout << "copy assignment" << std::endl;
    check(other.m_data);
    m_data = other.m_data;
    return *this;
  }
  Required& operator=(Required&& other) {
    std::cout << "move assignment" << std::endl;
    check(other.m_data);
    m_data = std::move(other.m_data);
    return *this;
  }

  T& operator*() noexcept { return get(); }
  const T& operator*() const noexcept { return get(); }
  T* operator->() noexcept { return &get(); }
  const T* operator->() const noexcept { return &get(); }

  constexpr explicit operator bool() const noexcept { return has(); }

  constexpr bool has() const noexcept { return m_data.has_value(); }

  constexpr T& get() & { return m_data.value(); }
  constexpr const T& get() const& { return m_data.value(); }
  constexpr T get() && { return m_data.value(); }

  template <class... Args>
  constexpr T& emplace(Args&&... args) {
    return m_data.emplace(std::forward<Args>(args)...);
  }

  template <class U, class... Args>
  constexpr T& emplace(std::initializer_list<U> ilist, Args&&... args) {
    return m_data.emplace(ilist, std::forward<Args>(args)...);
  }

 private:
  std::optional<T> m_data;

  static void check(const std::optional<T>& data) {
    if (!data) {
      throw std::runtime_error("data required");
    }
  }
};

}  // namespace ActsExamples

// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include <any>
#include <concepts>
#include <memory>
#include <stdexcept>
#include <type_traits>
#include <typeinfo>
#include <utility>

namespace detray {

/// Host-side owning type erasure for movable objects, including device buffers.
/// Const access is shallow, like unique_ptr: the owned object remains mutable.
class move_only_any {
 public:
  move_only_any() = default;
  move_only_any(const move_only_any&) = delete;
  move_only_any& operator=(const move_only_any&) = delete;
  move_only_any(move_only_any&&) noexcept = default;
  move_only_any& operator=(move_only_any&&) noexcept = default;

  template <typename T>
    requires(!std::same_as<std::remove_cvref_t<T>, move_only_any>)
  explicit move_only_any(T&& obj) {
    set(std::forward<T>(obj));
  }

  template <typename T>
    requires(!std::same_as<std::remove_cvref_t<T>, move_only_any>)
  void set(T&& obj) {
    using value_t = std::remove_cvref_t<T>;
    // Construct first: failed allocation/construction leaves the old value
    // intact. Typed new/delete also support over-aligned objects.
    storage_t replacement(new value_t(std::forward<T>(obj)),
                          [](void* ptr) { delete static_cast<value_t*>(ptr); });
    m_obj = std::move(replacement);
    m_type = &typeid(value_t);
  }

  void reset() noexcept { m_obj.reset(); }
  bool has_value() const noexcept { return static_cast<bool>(m_obj); }

  template <typename T>
  bool is() const noexcept {
    return has_value() && *m_type == typeid(T);
  }

  const std::type_info& type() const {
    if (!has_value()) {
      throw std::logic_error("Type requested from empty detray::move_only_any");
    }
    return *m_type;
  }

  template <typename T>
  T& as() const {
    if (!has_value()) {
      throw std::logic_error(
          "Value requested from empty detray::move_only_any");
    }
    if (!is<T>()) {
      throw std::bad_any_cast();
    }
    return *static_cast<T*>(m_obj.get());
  }

 private:
  using storage_t = std::unique_ptr<void, void (*)(void*)>;
  storage_t m_obj{nullptr, [](void*) {}};
  const std::type_info* m_type = nullptr;
};

}  // namespace detray

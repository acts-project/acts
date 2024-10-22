// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <array>
#include <cassert>
#include <cstddef>
#include <memory>

namespace Acts {

namespace detail {
/// Small opaque cache type which uses small buffer optimization
class SmallObjectCache {
 public:
  template <typename T, typename... Args>
  static SmallObjectCache make(Args&&... args) {
    SmallObjectCache cache{};

    static_assert(std::is_same_v<T, std::decay_t<T>>,
                  "Please pass the raw type, no const or ref");
    static_assert(sizeof(T) <= sizeof(cache.m_data),
                  "Passed type is too large");
    static_assert(
        std::is_move_assignable_v<T> && std::is_move_constructible_v<T>,
        "Type needs to be move assignable and move constructible");

    /*T* ptr =*/new (cache.m_data.data()) T(std::forward<Args>(args)...);
    static const Handler<T> static_handler{};
    cache.m_handler = &static_handler;

    return cache;
  }

  template <typename T>
  T& get() {
    static_assert(std::is_same_v<T, std::decay_t<T>>,
                  "Please pass the raw type, no const or ref");
    static_assert(sizeof(T) <= sizeof(m_data), "Passed type is too large");
    return *reinterpret_cast<T*>(m_data.data());
  }

  template <typename T>
  const T& get() const {
    static_assert(std::is_same_v<T, std::decay_t<T>>,
                  "Please pass the raw type, no const or ref");
    static_assert(sizeof(T) <= sizeof(m_data), "Passed type is too large");
    return *reinterpret_cast<const T*>(m_data.data());
  }

  ~SmallObjectCache() {
    assert(m_handler && "Handler cannot be dead");
    m_handler->destroy(m_data.data());
  }

  SmallObjectCache(SmallObjectCache&& other) {
    m_handler = other.m_handler;
    assert(m_handler && "Handler is null");
    m_handler->moveConstruct(other.m_data.data(), m_data.data());
  }

  SmallObjectCache& operator=(SmallObjectCache&& other) {
    m_handler = other.m_handler;
    assert(m_handler && "Handler is null");
    m_handler->move(other.m_data.data(), m_data.data());
    return *this;
  }

 private:
  SmallObjectCache() = default;

  struct HandlerBase {
    virtual void destroy(void* ptr) const = 0;
    virtual void moveConstruct(void* from, void* to) const = 0;
    virtual void move(void* from, void* to) const = 0;
    virtual ~HandlerBase() = default;
  };

  template <typename T>
  struct Handler final : public HandlerBase {
    void destroy(void* ptr) const override {
      assert(ptr != nullptr && "Address to destroy is nullptr");
      T* obj = static_cast<T*>(ptr);
      obj->~T();
    }

    void moveConstruct(void* from, void* to) const override {
      assert(from != nullptr && "Source is null");
      assert(to != nullptr && "Target is null");
      T* fromValue = static_cast<T*>(from);
      /*T* ptr =*/new (to) T(std::move(*fromValue));
    }

    void move(void* from, void* to) const override {
      assert(from != nullptr && "Source is null");
      assert(to != nullptr && "Target is null");

      T* fromValue = static_cast<T*>(from);
      T* toValue = static_cast<T*>(to);

      (*toValue) = std::move(*fromValue);
    }
  };

  alignas(std::max_align_t) std::array<char, 512> m_data{};
  const HandlerBase* m_handler{nullptr};
};

}  // namespace detail
}  // namespace Acts

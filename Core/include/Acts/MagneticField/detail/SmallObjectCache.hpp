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
template <size_t SIZE>
class SmallObjectCacheBase {
 public:
  template <typename T, typename... Args>
  static SmallObjectCacheBase make(Args&&... args) {
    SmallObjectCacheBase cache{};

    static_assert(std::is_same_v<T, std::decay_t<T>>,
                  "Please pass the raw type, no const or ref");
    static_assert(sizeof(T) <= SIZE, "Passed type is too large");
    static_assert(
        std::is_move_assignable_v<T> && std::is_move_constructible_v<T>,
        "Type needs to be move assignable and move constructible");

    /*T* ptr =*/new (cache.m_data.data()) T(std::forward<Args>(args)...);
    cache.m_handler = &makeHandler<T>();

    return cache;
  }

  template <typename T>
  explicit SmallObjectCacheBase(T&& value) {
    using U = std::decay_t<T>;
    static_assert(sizeof(T) <= SIZE, "Passed type is too large");
    static_assert(
        std::is_move_assignable_v<U> && std::is_move_constructible_v<U>,
        "Type needs to be move assignable and move constructible");

    // move construct into data block by move construct with placement new
    /*T* ptr =*/new (m_data.data()) U(std::move(value));
    m_handler = &makeHandler<U>();
  }

  template <typename T>
  T& get() {
    static_assert(std::is_same_v<T, std::decay_t<T>>,
                  "Please pass the raw type, no const or ref");
    static_assert(sizeof(T) <= SIZE, "Passed type is too large");
    return *reinterpret_cast<T*>(m_data.data());
  }

  template <typename T>
  const T& get() const {
    static_assert(std::is_same_v<T, std::decay_t<T>>,
                  "Please pass the raw type, no const or ref");
    static_assert(sizeof(T) <= SIZE, "Passed type is too large");
    return *reinterpret_cast<const T*>(m_data.data());
  }

  ~SmallObjectCacheBase() {
    assert(m_handler && "Handler cannot be dead");
    if (m_handler != nullptr) {
      m_handler->destroy(m_data.data());
    }
  }

  SmallObjectCacheBase(const SmallObjectCacheBase& other) {
    if (m_handler == nullptr) {
      m_handler = other.m_handler;
      assert(m_handler != nullptr && "Handler is null");
      m_handler->copyConstruct(other.m_data.data(), m_data.data());
    } else {
      assert(m_handler == other.m_handler);
      m_handler->copy(other.m_data.data(), m_data.data());
    }
  }

  SmallObjectCacheBase& operator=(const SmallObjectCacheBase& other) {
    if (m_handler == nullptr) {
      m_handler = other.m_handler;
      assert(m_handler && "Handler is null");
      m_handler->copyConstruct(other.m_data.data(), m_data.data());
    } else {
      assert(m_handler == other.m_handler);
      m_handler->copy(other.m_data.data(), m_data.data());
    }
    return *this;
  }

  SmallObjectCacheBase(SmallObjectCacheBase&& other) {
    if (m_handler == nullptr) {
      m_handler = other.m_handler;
      assert(m_handler != nullptr && "Handler is null");
      m_handler->moveConstruct(other.m_data.data(), m_data.data());
    } else {
      assert(m_handler == other.m_handler);
      m_handler->move(other.m_data.data(), m_data.data());
    }
  }

  SmallObjectCacheBase& operator=(SmallObjectCacheBase&& other) {
    if (m_handler == nullptr) {
      m_handler = other.m_handler;
      assert(m_handler && "Handler is null");
      m_handler->moveConstruct(other.m_data.data(), m_data.data());
    } else {
      assert(m_handler == other.m_handler);
      m_handler->move(other.m_data.data(), m_data.data());
    }
    return *this;
  }

  SmallObjectCacheBase() = default;

 private:
  struct Handler {
    void (*destroy)(void* ptr);
    void (*moveConstruct)(void* from, void* to);
    void (*move)(void* from, void* to);
    void (*copyConstruct)(const void* from, void* to);
    void (*copy)(const void* from, void* to);

    // virtual void destroy(void* ptr) const = 0;
    // virtual void moveConstruct(void* from, void* to) const = 0;
    // virtual void move(void* from, void* to) const = 0;
    // virtual void copyConstruct(const void* from, void* to) const = 0;
    // virtual void copy(const void* from, void* to) const = 0;
    // virtual ~HandlerBase() = default;
  };

  template <typename T>
  static Handler& makeHandler() {
    static Handler static_handler = []() {
      Handler h;
      h.destroy = &destroy<T>;
      h.moveConstruct = &moveConstruct<T>;
      h.move = &move<T>;
      h.copyConstruct = &copyConstruct<T>;
      h.copy = &copy<T>;
      return h;
    }();
    return static_handler;
  }

  template <typename T>
  static void destroy(void* ptr) {
    assert(ptr != nullptr && "Address to destroy is nullptr");
    T* obj = static_cast<T*>(ptr);
    obj->~T();
  }

  template <typename T>
  static void moveConstruct(void* from, void* to) {
    assert(from != nullptr && "Source is null");
    assert(to != nullptr && "Target is null");
    T* _from = static_cast<T*>(from);
    /*T* ptr =*/new (to) T(std::move(*_from));
  }

  template <typename T>
  static void move(void* from, void* to) {
    assert(from != nullptr && "Source is null");
    assert(to != nullptr && "Target is null");

    T* _from = static_cast<T*>(from);
    T* _to = static_cast<T*>(to);

    (*_to) = std::move(*_from);
  }

  template <typename T>
  static void copyConstruct(const void* from, void* to) {
    assert(from != nullptr && "Source is null");
    assert(to != nullptr && "Target is null");
    const T* _from = static_cast<const T*>(from);
    /*T* ptr =*/new (to) T(*_from);
  }

  template <typename T>
  static void copy(const void* from, void* to) {
    assert(from != nullptr && "Source is null");
    assert(to != nullptr && "Target is null");

    const T* _from = static_cast<const T*>(from);
    T* _to = static_cast<T*>(to);

    (*_to) = *_from;
  }

  alignas(std::max_align_t) std::array<char, SIZE> m_data{};
  Handler* m_handler{nullptr};
};

using SmallObjectCache = SmallObjectCacheBase<512>;

}  // namespace detail
}  // namespace Acts

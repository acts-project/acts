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

/// Small opaque cache type which uses small buffer optimization
template <size_t SIZE>
class AnyBase {
  static_assert(SIZE > sizeof(void*), "Size cannot be smaller than a pointer");

 public:
  template <typename T, typename... Args>
  static AnyBase make(Args&&... args) {
    AnyBase cache{};

    using U = std::decay_t<T>;
    static_assert(std::is_same_v<T, std::decay_t<T>>,
                  "Please pass the raw type, no const or ref");
    static_assert(sizeof(U) <= sizeof(m_data), "Passed type is too large");
    static_assert(
        std::is_move_assignable_v<U> && std::is_move_constructible_v<U>,
        "Type needs to be move assignable and move constructible");
    static_assert(
        std::is_copy_assignable_v<U> && std::is_copy_constructible_v<U>,
        "Type needs to be copy assignable and copy constructible");

    /*T* ptr =*/new (cache.m_data.data()) U(std::forward<Args>(args)...);
    cache.m_handler = makeHandler<U>();

    return cache;
  }

  template <typename T>
  explicit AnyBase(T&& value) {
    using U = std::decay_t<T>;
    static_assert(sizeof(U) <= sizeof(m_data), "Passed type is too large");
    static_assert(
        std::is_move_assignable_v<U> && std::is_move_constructible_v<U>,
        "Type needs to be move assignable and move constructible");
    static_assert(
        std::is_copy_assignable_v<U> && std::is_copy_constructible_v<U>,
        "Type needs to be copy assignable and copy constructible");

    // move construct into data block by move construct with placement new
    /*T* ptr =*/new (m_data.data()) U(std::move(value));
    m_handler = makeHandler<U>();
  }

  template <typename T>
  T& get() {
    static_assert(std::is_same_v<T, std::decay_t<T>>,
                  "Please pass the raw type, no const or ref");
    static_assert(sizeof(T) <= sizeof(m_data), "Passed type is too large");
    assert(makeHandler<T>() == m_handler && "Bad type access");
    return *reinterpret_cast<T*>(m_data.data());
  }

  template <typename T>
  const T& get() const {
    static_assert(std::is_same_v<T, std::decay_t<T>>,
                  "Please pass the raw type, no const or ref");
    static_assert(sizeof(T) <= sizeof(m_data), "Passed type is too large");
    assert(makeHandler<T>() == m_handler && "Bad type access");
    return *reinterpret_cast<const T*>(m_data.data());
  }

  ~AnyBase() {
    assert(m_handler && "Handler cannot be dead");
    destroy();
  }

  AnyBase(const AnyBase& other) {
    if (m_handler == nullptr) {  // this object is empty
      m_handler = other.m_handler;
      assert(m_handler != nullptr && "Handler is null");
      copyConstruct(other);
    } else {
      // @TODO: Support assigning between different types
      assert(m_handler == other.m_handler);
      copy(other);
    }
  }

  AnyBase& operator=(const AnyBase& other) {
    if (m_handler == nullptr) {  // this object is empty
      m_handler = other.m_handler;
      assert(m_handler && "Handler is null");
      copyConstruct(other);
    } else {
      // @TODO: Support assigning between different types
      assert(m_handler == other.m_handler);
      copy(other);
    }
    return *this;
  }

  AnyBase(AnyBase&& other) {
    if (m_handler == nullptr) {  // this object is empty
      m_handler = other.m_handler;
      assert(m_handler != nullptr && "Handler is null");
      moveConstruct(std::move(other));
    } else {
      // @TODO: Support assigning between different types
      assert(m_handler == other.m_handler);
      move(std::move(other));
    }
  }

  AnyBase& operator=(AnyBase&& other) {
    if (m_handler == nullptr) {  // this object is empty
      m_handler = other.m_handler;
      assert(m_handler && "Handler is null");
      moveConstruct(std::move(other));
    } else {
      // @TODO: Support assigning between different types
      assert(m_handler == other.m_handler);
      move(std::move(other));
    }
    return *this;
  }

  AnyBase() = default;

 private:
  struct Handler {
    void (*destroy)(void* ptr) = nullptr;
    void (*moveConstruct)(void* from, void* to) = nullptr;
    void (*move)(void* from, void* to) = nullptr;
    void (*copyConstruct)(const void* from, void* to) = nullptr;
    void (*copy)(const void* from, void* to) = nullptr;
  };

  template <typename T>
  static Handler* makeHandler() {
    static Handler static_handler = []() {
      Handler h;
      if constexpr (!std::is_trivially_destructible_v<T>) {
        h.destroy = &destroyImpl<T>;
      }
      if constexpr (!std::is_trivially_move_constructible_v<T>) {
        h.moveConstruct = &moveConstructImpl<T>;
      }
      if constexpr (!std::is_trivially_move_assignable_v<T>) {
        h.move = &moveImpl<T>;
      }
      if constexpr (!std::is_trivially_copy_constructible_v<T>) {
        h.copyConstruct = &copyConstructImpl<T>;
      }
      if constexpr (!std::is_trivially_copy_assignable_v<T>) {
        h.copy = &copyImpl<T>;
      }
      return h;
    }();
    return &static_handler;
  }

  void destroy() {
    if (m_handler != nullptr && m_handler->destroy != nullptr) {
      m_handler->destroy(m_data.data());
    }
  }

  void moveConstruct(AnyBase&& from) {
    if (m_handler == nullptr) {
      return;
    }
    if (m_handler->moveConstruct == nullptr) {
      // trivially move constructible
      m_data = std::move(from.m_data);
    } else {
      m_handler->moveConstruct(from.m_data.data(), m_data.data());
    }
  }

  void move(AnyBase&& from) {
    if (m_handler == nullptr) {
      return;
    }
    if (m_handler->move == nullptr) {
      // trivially move constructible -> trivially movable
      m_data = std::move(from.m_data);
    } else {
      m_handler->move(from.m_data.data(), m_data.data());
    }
  }

  void copyConstruct(const AnyBase& from) {
    if (m_handler == nullptr) {
      return;
    }
    if (m_handler->copyConstruct == nullptr) {
      // trivially copyable
      m_data = from.m_data;
    } else {
      m_handler->copyConstruct(from.m_data.data(), m_data.data());
    }
  }

  void copy(const AnyBase& from) {
    if (m_handler == nullptr) {
      return;
    }
    if (m_handler->copy == nullptr) {
      // trivially copyable
      m_data = from.m_data;
    } else {
      m_handler->copy(from.m_data.data(), m_data.data());
    }
  }

  template <typename T>
  static void destroyImpl(void* ptr) {
    assert(ptr != nullptr && "Address to destroy is nullptr");
    T* obj = static_cast<T*>(ptr);
    obj->~T();
  }

  template <typename T>
  static void moveConstructImpl(void* from, void* to) {
    assert(from != nullptr && "Source is null");
    assert(to != nullptr && "Target is null");
    T* _from = static_cast<T*>(from);
    /*T* ptr =*/new (to) T(std::move(*_from));
  }

  template <typename T>
  static void moveImpl(void* from, void* to) {
    assert(from != nullptr && "Source is null");
    assert(to != nullptr && "Target is null");

    T* _from = static_cast<T*>(from);
    T* _to = static_cast<T*>(to);

    (*_to) = std::move(*_from);
  }

  template <typename T>
  static void copyConstructImpl(const void* from, void* to) {
    assert(from != nullptr && "Source is null");
    assert(to != nullptr && "Target is null");
    const T* _from = static_cast<const T*>(from);
    /*T* ptr =*/new (to) T(*_from);
  }

  template <typename T>
  static void copyImpl(const void* from, void* to) {
    assert(from != nullptr && "Source is null");
    assert(to != nullptr && "Target is null");

    const T* _from = static_cast<const T*>(from);
    T* _to = static_cast<T*>(to);

    (*_to) = *_from;
  }

  alignas(std::max_align_t) std::array<char, SIZE> m_data{};
  Handler* m_handler{nullptr};
};

using Any = AnyBase<8>;

}  // namespace Acts

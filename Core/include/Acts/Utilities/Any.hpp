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
#include <utility>

// #define _ENABLE_ACTS_ANY_VERBOSE
// #define _ENABLE_ACTS_ANY_DEBUG
// #define _ACTS_ANY_TRACK_ALLOCATIONS

#if defined(_ACTS_ANY_TRACK_ALLOCATIONS)
#include <iostream>
#include <mutex>
#include <set>
#include <typeindex>
#include <typeinfo>
#endif

#if defined(_ENABLE_ACTS_ANY_VERBOSE) || defined(_ENABLE_ACTS_ANY_DEBUG)
#include <iomanip>
#include <iostream>
#endif

#if defined(_ENABLE_ACTS_ANY_DEBUG)
#define _ACTS_ANY_DEBUG(x) std::cout << x << std::endl;
#else
#define _ACTS_ANY_DEBUG(x)
#endif

#if defined(_ENABLE_ACTS_ANY_VERBOSE)
#define _ACTS_ANY_VERBOSE(x) std::cout << x << std::endl;
#define _ACTS_ANY_VERBOSE_BUFFER(s, b) \
  do {                                 \
    std::cout << "" << s << ": 0x";    \
    for (char c : b) {                 \
      std::cout << std::hex << (int)c; \
    }                                  \
    std::cout << std::endl;            \
  } while (0)
#else
#define _ACTS_ANY_VERBOSE(x)
#define _ACTS_ANY_VERBOSE_BUFFER(s, b)
#endif

namespace Acts {

#if defined(_ACTS_ANY_TRACK_ALLOCATIONS)
static std::mutex _s_any_mutex;
static std::set<std::pair<std::type_index, void*>> _s_any_allocations;

struct _AnyAllocationReporter {
  ~_AnyAllocationReporter() {
    std::lock_guard guard{_s_any_mutex};

    if (!_s_any_allocations.empty()) {
      std::cout << "Not all allocations have been released" << std::endl;
      std::abort();
    }
  }
};
static _AnyAllocationReporter s_reporter;
#endif

class AnyBaseAll {};

/// Small opaque cache type which uses small buffer optimization
template <size_t SIZE>
class AnyBase : public AnyBaseAll {
  static_assert(sizeof(void*) <= SIZE, "Size is too small for a pointer");

 public:
  template <typename T, typename... Args>
  static AnyBase make(Args&&... args) {
    AnyBase cache{};

    using U = std::decay_t<T>;
    static_assert(std::is_same_v<T, std::decay_t<T>>,
                  "Please pass the raw type, no const or ref");
    static_assert(
        std::is_move_assignable_v<U> && std::is_move_constructible_v<U>,
        "Type needs to be move assignable and move constructible");
    static_assert(
        std::is_copy_assignable_v<U> && std::is_copy_constructible_v<U>,
        "Type needs to be copy assignable and copy constructible");

    cache.m_handler = makeHandler<U>();
    if constexpr (sizeof(U) <= SIZE) {
      // construct into local buffer
      /*U* ptr =*/new (cache.m_data.data()) U(std::forward<Args>(args)...);
    } else {
      // too large, heap allocate
      U* heap = new U(std::forward<Args>(args)...);
#if defined(_ACTS_ANY_TRACK_ALLOCATIONS)
      {
        std::lock_guard guard{_s_any_mutex};
        _s_any_allocations.emplace(std::type_index(typeid(T)), heap);
      }
#endif
      *reinterpret_cast<U**>(cache.m_data.data()) = heap;
    }

    return cache;
  }

  template <typename T, typename... Args>
  explicit AnyBase(std::in_place_type_t<T> /*unused*/, Args&&... args) {
    using U = std::decay_t<T>;
    static_assert(std::is_same_v<T, std::decay_t<T>>,
                  "Please pass the raw type, no const or ref");
    static_assert(
        std::is_move_assignable_v<U> && std::is_move_constructible_v<U>,
        "Type needs to be move assignable and move constructible");
    static_assert(
        std::is_copy_assignable_v<U> && std::is_copy_constructible_v<U>,
        "Type needs to be copy assignable and copy constructible");

    m_handler = makeHandler<U>();
    if constexpr (sizeof(U) <= SIZE) {
      // construct into local buffer
      /*U* ptr =*/new (m_data.data()) U(std::forward<Args>(args)...);
    } else {
      // too large, heap allocate
      U* heap = new U(std::forward<Args>(args)...);
#if defined(_ACTS_ANY_TRACK_ALLOCATIONS)
      {
        std::lock_guard guard{_s_any_mutex};
        _s_any_allocations.emplace(std::type_index(typeid(T)), heap);
      }
#endif
      *reinterpret_cast<U**>(m_data.data()) = heap;
    }
  }

  AnyBase() {
    m_data.fill(0);
#if defined(_ENABLE_ACTS_ANY_VERBOSE)
    _ACTS_ANY_VERBOSE("Default construct this=" << this);
#endif
  };

  template <typename T, typename = std::enable_if_t<
                            !std::is_same_v<std::decay_t<T>, AnyBase<SIZE>>>>
  explicit AnyBase(T&& value) {
    using U = std::decay_t<T>;
    static_assert(
        std::is_move_assignable_v<U> && std::is_move_constructible_v<U>,
        "Type needs to be move assignable and move constructible");
    static_assert(
        std::is_copy_assignable_v<U> && std::is_copy_constructible_v<U>,
        "Type needs to be copy assignable and copy constructible");

    m_handler = makeHandler<U>();

    if constexpr (sizeof(U) <= SIZE) {
      // construct into local buffer
      /*U* ptr =*/new (m_data.data()) U(std::move(value));
      _ACTS_ANY_VERBOSE(
          "Construct local (this=" << this << ") at: " << (void*)m_data.data());
    } else {
      // too large, heap allocate
      U* heap = new U(std::move(value));
      _ACTS_ANY_VERBOSE("Construct heap (this=" << this << ") at: " << heap);
      _ACTS_ANY_VERBOSE_BUFFER("-> buffer before", m_data);
#if defined(_ACTS_ANY_TRACK_ALLOCATIONS)
      {
        std::lock_guard guard{_s_any_mutex};
        _s_any_allocations.emplace(std::type_index(typeid(T)), heap);
      }
#endif
      *reinterpret_cast<U**>(m_data.data()) = heap;
      _ACTS_ANY_VERBOSE_BUFFER("-> buffer after", m_data);
    }
  }

  // template <typename T>
  // explicit AnyBase(const T& value) {
  // using U = std::decay_t<T>;
  // static_assert(
  // std::is_move_assignable_v<U> && std::is_move_constructible_v<U>,
  // "Type needs to be move assignable and move constructible");
  // static_assert(
  // std::is_copy_assignable_v<U> && std::is_copy_constructible_v<U>,
  // "Type needs to be copy assignable and copy constructible");

  // m_handler = makeHandler<U>();

  // if constexpr (sizeof(U) <= SIZE) {
  // // construct into local buffer
  // [>U* ptr =<]new (m_data.data()) U(value);
  // _ACTS_ANY_VERBOSE(
  // "Construct local (this=" << this << ") at: " << (void*)m_data.data());
  // } else {
  // // too large, heap allocate
  // U* heap = new U(value);
  // _ACTS_ANY_VERBOSE("Construct heap (this=" << this << ") at: " << heap);
  // _ACTS_ANY_VERBOSE_BUFFER("-> buffer before", m_data);
  // #if defined(_ACTS_ANY_TRACK_ALLOCATIONS)
  // {
  // std::lock_guard guard{_s_any_mutex};
  // _s_any_allocations.emplace(std::type_index(typeid(T)), heap);
  // }
  // #endif
  // *reinterpret_cast<U**>(m_data.data()) = heap;
  // _ACTS_ANY_VERBOSE_BUFFER("-> buffer after", m_data);
  // }
  // }

  template <typename T>
  T& as() {
    static_assert(std::is_same_v<T, std::decay_t<T>>,
                  "Please pass the raw type, no const or ref");
    assert(makeHandler<T>() == m_handler && "Bad type access");
    if (m_handler->typeSize <= SIZE) {
      T* ptr = reinterpret_cast<T*>(m_data.data());
      _ACTS_ANY_VERBOSE("As local: " << ptr);
      return *ptr;
    } else {
      // is heap allocated
      void* ptr = *reinterpret_cast<void**>(m_data.data());
      _ACTS_ANY_VERBOSE("As heap: " << ptr);
      return *reinterpret_cast<T*>(ptr);
    }
  }

  template <typename T>
  const T& as() const {
    static_assert(std::is_same_v<T, std::decay_t<T>>,
                  "Please pass the raw type, no const or ref");
    assert(makeHandler<T>() == m_handler && "Bad type access");
    if (m_handler->typeSize <= SIZE) {
      const T* ptr = reinterpret_cast<const T*>(m_data.data());
      _ACTS_ANY_VERBOSE("As local: " << ptr);
      return *ptr;
    } else {
      // is heap allocated
      const void* ptr = *reinterpret_cast<void* const*>(m_data.data());
      _ACTS_ANY_VERBOSE("As heap: " << ptr);
      return *reinterpret_cast<const T*>(ptr);
    }
  }

  ~AnyBase() { destroy(); }

  AnyBase(const AnyBase& other) {
    if (m_handler == nullptr) {  // this object is empty
      m_handler = other.m_handler;
      assert(m_handler != nullptr && "Handler is null");
      // @FIXME: Need to allocate
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
      // @FIXME: Need to allocate
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

  operator bool() const { return m_handler != nullptr; }

 private:
  struct Handler {
    void (*destroy)(void* ptr) = nullptr;
    void (*moveConstruct)(void* from, void* to) = nullptr;
    void (*move)(void* from, void* to) = nullptr;
    void* (*copyConstruct)(const void* from, void* to) = nullptr;
    void (*copy)(const void* from, void* to) = nullptr;
    std::size_t typeSize{0};
  };

  template <typename T>
  static Handler* makeHandler() {
    static_assert(!std::is_same_v<T, AnyBase<SIZE>>, "Cannot wrap any in any");
    static Handler static_handler = []() {
      Handler h;
      if constexpr (!std::is_trivially_destructible_v<T> || sizeof(T) > SIZE) {
        h.destroy = &destroyImpl<T>;
      }
      if constexpr (!std::is_trivially_move_constructible_v<T> ||
                    sizeof(T) > SIZE) {
        h.moveConstruct = &moveConstructImpl<T>;
      }
      if constexpr (!std::is_trivially_move_assignable_v<T> ||
                    sizeof(T) > SIZE) {
        h.move = &moveImpl<T>;
      }
      if constexpr (!std::is_trivially_copy_constructible_v<T> ||
                    sizeof(T) > SIZE) {
        h.copyConstruct = &copyConstructImpl<T>;
      }
      if constexpr (!std::is_trivially_copy_assignable_v<T> ||
                    sizeof(T) > SIZE) {
        h.copy = &copyImpl<T>;
      }
      h.typeSize = sizeof(T);

      _ACTS_ANY_DEBUG("Type: " << typeid(T).name());
      _ACTS_ANY_DEBUG(" -> destroy: " << h.destroy);
      _ACTS_ANY_DEBUG(" -> moveConstruct: " << h.moveConstruct);
      _ACTS_ANY_DEBUG(" -> move: " << h.move);
      _ACTS_ANY_DEBUG(" -> copyConstruct: " << h.copyConstruct);
      _ACTS_ANY_DEBUG(" -> copy: " << h.copy);
      _ACTS_ANY_DEBUG(" -> typeSize: " << h.typeSize << " heap: "
                                       << (h.typeSize > SIZE ? "yes" : "no"));

      return h;
    }();
    return &static_handler;
  }

  void destroy() {
    _ACTS_ANY_VERBOSE("Destructor this=" << this << " handler: " << m_handler);
    if (m_handler != nullptr && m_handler->destroy != nullptr) {
      void* ptr = m_data.data();
      if (m_handler->typeSize > SIZE) {
        // stored on heap: interpret buffer as pointer
        _ACTS_ANY_VERBOSE("Pre-destroy data: " << (void*)m_data.data());

        _ACTS_ANY_VERBOSE_BUFFER("-> buffer pre-destroy", m_data);

        ptr = *reinterpret_cast<void**>(m_data.data());
        _ACTS_ANY_VERBOSE("Pre-destroy: " << ptr);
      }
      m_handler->destroy(ptr);
      m_handler = nullptr;
    }
  }

  void moveConstruct(AnyBase&& fromAny) {
    if (m_handler == nullptr) {
      return;
    }

    void* to = m_data.data();
    void* from = fromAny.m_data.data();
    if (m_handler->typeSize > SIZE) {
      // stored on heap: just copy the pointer
      *reinterpret_cast<void**>(m_data.data()) =
          *reinterpret_cast<void**>(fromAny.m_data.data());
      // do not delete in moved-from any
      fromAny.m_handler = nullptr;
      return;
    }

    if (m_handler->moveConstruct == nullptr) {
      // trivially move constructible
      m_data = std::move(fromAny.m_data);
    } else {
      m_handler->moveConstruct(from, to);
    }
  }

  void move(AnyBase&& fromAny) {
    if (m_handler == nullptr) {
      return;
    }

    void* to = m_data.data();
    void* from = fromAny.m_data.data();
    if (m_handler->typeSize > SIZE) {
      // stored on heap: just copy the pointer
      *reinterpret_cast<void**>(m_data.data()) =
          *reinterpret_cast<void**>(fromAny.m_data.data());
      // do not delete in moved-from any
      fromAny.m_handler = nullptr;
      return;
    }

    if (m_handler->move == nullptr) {
      // trivially move constructible -> trivially movable
      m_data = std::move(fromAny.m_data);
    } else {
      m_handler->move(from, to);
    }
  }

  void copyConstruct(const AnyBase& fromAny) {
    if (m_handler == nullptr) {
      return;
    }

    void* to = m_data.data();
    const void* from = fromAny.m_data.data();
    if (m_handler->typeSize > SIZE) {
      // stored on heap: interpret buffer as pointer
      to = *reinterpret_cast<void**>(m_data.data());
      from = *reinterpret_cast<void* const*>(fromAny.m_data.data());

      // if (to == nullptr) {
      // // want to copy but do not have allocation yet
      // to = new char[m_handler->typeSize];
      // *reinterpret_cast<void**>(m_data.data()) = to;
      // }
    }

    if (m_handler->copyConstruct == nullptr) {
      // trivially copyable
      m_data = fromAny.m_data;
    } else {
      void* copyAt = m_handler->copyConstruct(from, to);
      if (to == nullptr) {
        assert(copyAt != nullptr);
        // copy allocated, store pointer
        *reinterpret_cast<void**>(m_data.data()) = copyAt;
      }
    }
  }

  void copy(const AnyBase& fromAny) {
    if (m_handler == nullptr) {
      return;
    }

    void* to = m_data.data();
    const void* from = fromAny.m_data.data();
    if (m_handler->typeSize > SIZE) {
      // stored on heap: interpret buffer as pointer
      to = *reinterpret_cast<void**>(m_data.data());
      from = *reinterpret_cast<void* const*>(fromAny.m_data.data());
    }

    if (m_handler->copy == nullptr) {
      // trivially copyable
      m_data = fromAny.m_data;
    } else {
      m_handler->copy(from, to);
    }
  }

  template <typename T>
  static void destroyImpl(void* ptr) {
    assert(ptr != nullptr && "Address to destroy is nullptr");
    T* obj = static_cast<T*>(ptr);
    if constexpr (sizeof(T) <= SIZE) {
      // stored in place: just call the destructor
      _ACTS_ANY_VERBOSE("Destroy local at: " << ptr);
      obj->~T();
    } else {
      // stored on heap: delete
      _ACTS_ANY_VERBOSE("Delete heap at: " << obj);
#if defined(_ACTS_ANY_TRACK_ALLOCATIONS)
      {
        std::lock_guard guard{_s_any_mutex};
        auto it =
            _s_any_allocations.find(std::pair{std::type_index(typeid(T)), obj});
        if (it == _s_any_allocations.end()) {
          throw std::runtime_error{
              "Trying to deallocate heap address that we didn't allocate"};
        }
        _s_any_allocations.erase(it);
      }
#endif
      delete obj;
    }
  }

  template <typename T>
  static void moveConstructImpl(void* from, void* to) {
    _ACTS_ANY_VERBOSE("move const: " << from << " -> " << to);
    assert(from != nullptr && "Source is null");
    assert(to != nullptr && "Target is null");
    T* _from = static_cast<T*>(from);
    /*T* ptr =*/new (to) T(std::move(*_from));
  }

  template <typename T>
  static void moveImpl(void* from, void* to) {
    _ACTS_ANY_VERBOSE("move: " << from << " -> " << to);
    assert(from != nullptr && "Source is null");
    assert(to != nullptr && "Target is null");

    T* _from = static_cast<T*>(from);
    T* _to = static_cast<T*>(to);

    (*_to) = std::move(*_from);
  }

  template <typename T>
  static void* copyConstructImpl(const void* from, void* to) {
    _ACTS_ANY_VERBOSE("copy const: " << from << " -> " << to);
    assert(from != nullptr && "Source is null");
    // assert(to != nullptr && "Target is null");
    const T* _from = static_cast<const T*>(from);
    if (to == nullptr) {
      assert(sizeof(T) > SIZE && "Received nullptr in local buffer case");
      to = new T(*_from);
#if defined(_ACTS_ANY_TRACK_ALLOCATIONS)
      {
        std::lock_guard guard{_s_any_mutex};
        _s_any_allocations.emplace(std::type_index(typeid(T)), to);
      }
#endif

    } else {
      assert(sizeof(T) <= SIZE && "Received non-nullptr in heap case");
      /*T* ptr =*/new (to) T(*_from);
    }
    return to;
  }

  template <typename T>
  static void copyImpl(const void* from, void* to) {
    _ACTS_ANY_VERBOSE("copy: " << from << " -> " << to);
    assert(from != nullptr && "Source is null");
    assert(to != nullptr && "Target is null");

    const T* _from = static_cast<const T*>(from);
    T* _to = static_cast<T*>(to);

    (*_to) = *_from;
  }

  alignas(std::max_align_t) std::array<char, SIZE> m_data{};
  const Handler* m_handler{nullptr};
};

#if defined(_ACTS_ANY_TRACK_ALLOCATIONS)
#endif

using Any = AnyBase<8>;

#undef _ACTS_ANY_VERBOSE
#undef _ACTS_ANY_VERBOSE_BUFFER
#undef _ENABLE_ACTS_ANY_VERBOSE

}  // namespace Acts

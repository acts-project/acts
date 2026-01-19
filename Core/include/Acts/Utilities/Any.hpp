// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Utilities/HashedString.hpp"

#include <any>
#include <array>
#include <cassert>
#include <cstddef>
#include <utility>

// #define _ACTS_ANY_ENABLE_VERBOSE
// #define _ACTS_ANY_ENABLE_DEBUG
// #define _ACTS_ANY_ENABLE_TRACK_ALLOCATIONS

#if defined(_ACTS_ANY_ENABLE_TRACK_ALLOCATIONS)
#include <iostream>
#include <mutex>
#include <set>
#include <typeindex>
#include <typeinfo>
#endif

#if defined(_ACTS_ANY_ENABLE_VERBOSE) || defined(_ACTS_ANY_ENABLE_DEBUG)
#include <iomanip>
#include <iostream>
#endif

#if defined(_ACTS_ANY_ENABLE_DEBUG)
#define _ACTS_ANY_DEBUG(x) std::cout << x << std::endl;
#else
#define _ACTS_ANY_DEBUG(x)
#endif

#if defined(_ACTS_ANY_ENABLE_VERBOSE)
#define _ACTS_ANY_VERBOSE(x) std::cout << x << std::endl;
#define _ACTS_ANY_VERBOSE_BUFFER(s, b)              \
  do {                                              \
    std::cout << "" << s << ": 0x";                 \
    for (char c : b) {                              \
      std::cout << std::hex << static_cast<int>(c); \
    }                                               \
    std::cout << std::endl;                         \
  } while (0)
#else
#define _ACTS_ANY_VERBOSE(x)
#define _ACTS_ANY_VERBOSE_BUFFER(s, b)
#endif

namespace Acts {

#if defined(_ACTS_ANY_ENABLE_TRACK_ALLOCATIONS)
static std::mutex _s_any_mutex;
static std::set<std::pair<std::type_index, void*>> _s_any_allocations;

#define _ACTS_ANY_TRACK_ALLOCATION(T, heap)                                  \
  do {                                                                       \
    std::lock_guard guard{_s_any_mutex};                                     \
    _s_any_allocations.emplace(std::type_index(typeid(T)), heap);            \
    _ACTS_ANY_DEBUG("Allocate type: " << typeid(T).name() << " at " << heap) \
  } while (0)

#define _ACTS_ANY_TRACK_DEALLOCATION(T, heap)                                 \
  do {                                                                        \
    std::lock_guard guard{_s_any_mutex};                                      \
    auto it =                                                                 \
        _s_any_allocations.find(std::pair{std::type_index(typeid(T)), heap}); \
    if (it == _s_any_allocations.end()) {                                     \
      throw std::runtime_error{                                               \
          "Trying to deallocate heap address that we didn't allocate"};       \
    }                                                                         \
    _s_any_allocations.erase(it);                                             \
  } while (0)

// Do not make member functions noexcept in the debug case
#define _ACTS_ANY_NOEXCEPT /*nothing*/

struct _AnyAllocationReporter {
  static void checkAllocations() {
    std::lock_guard guard{_s_any_mutex};

    if (!_s_any_allocations.empty()) {
      std::cout << "Not all allocations have been released" << std::endl;
      for (const auto& [idx, addr] : _s_any_allocations) {
        std::cout << "- " << idx.name() << ": " << addr << std::endl;
      }
      throw std::runtime_error{"AnyCheckAllocations failed"};
    }
  }

  ~_AnyAllocationReporter() noexcept { checkAllocations(); }
};
static _AnyAllocationReporter s_reporter;
#else
#define _ACTS_ANY_TRACK_ALLOCATION(T, heap) \
  do {                                      \
  } while (0)
#define _ACTS_ANY_TRACK_DEALLOCATION(T, heap) \
  do {                                        \
  } while (0)
#define _ACTS_ANY_NOEXCEPT noexcept
#endif

/// @addtogroup utilities
/// @{

/// Base class for all instances of @ref AnyBase regarfless of SBO size
class AnyBaseAll {};

/// Small opaque type-erased type with configurable small buffer optimization
///
/// @note
/// Type requirements:
/// - All stored types must be copy constructible and copy assignable.
/// - Types stored locally (`sizeof(T) <= SIZE`) must also be move constructible
///   and move assignable because local moves use move operations when not
///   trivially movable (trivial moves fall back to buffer copies).
/// - Types stored on the heap (`sizeof(T) > SIZE`) are moved by stealing the
///   pointer, so no move operations are required in that case.
///
/// @note
/// In summary:
/// - Local storage: values live inside the internal buffer; moves may invoke
///   move operations or buffer copies; copies use copy operations or buffer
///   copies when trivial.
/// - Heap storage: values are allocated on the heap; moves transfer ownership
///   of the pointer; copies allocate and copy-construct the pointee.
template <std::size_t SIZE>
class AnyBase : public AnyBaseAll {
  static_assert(sizeof(void*) <= SIZE, "Size is too small for a pointer");

 public:
  /// Construct with in-place type construction
  /// @tparam T Type to construct
  /// @tparam Args Constructor argument types
  /// @param args Arguments to forward to T's constructor
  template <typename T, typename... Args>
    requires(std::is_copy_assignable_v<std::decay_t<T>> &&
             std::is_copy_constructible_v<std::decay_t<T>> &&
             (sizeof(std::decay_t<T>) > SIZE ||
              (std::is_move_assignable_v<std::decay_t<T>> &&
               std::is_move_constructible_v<std::decay_t<T>>)))
  explicit AnyBase(std::in_place_type_t<T> /*unused*/, Args&&... args) {
    using U = std::decay_t<T>;
    m_handler = makeHandler<U>();
    constructValue<U>(std::forward<Args>(args)...);
  }

#if defined(_ACTS_ANY_ENABLE_VERBOSE)
  AnyBase() { _ACTS_ANY_VERBOSE("Default construct this=" << this); };
#else
  AnyBase() = default;
#endif

  /// Construct from any value type
  /// @tparam T Type of the value to store
  /// @param value Value to store in the Any
  template <typename T>
  explicit AnyBase(T&& value) _ACTS_ANY_NOEXCEPT
    requires(!std::same_as<std::decay_t<T>, AnyBase<SIZE>> &&
             std::is_copy_assignable_v<std::decay_t<T>> &&
             std::is_copy_constructible_v<std::decay_t<T>> &&
             (sizeof(std::decay_t<T>) > SIZE ||
              (std::is_move_assignable_v<std::decay_t<T>> &&
               std::is_move_constructible_v<std::decay_t<T>>)))
      : AnyBase{std::in_place_type<T>, std::forward<T>(value)} {}

  /// Construct a new value in place, destroying any existing value
  /// @tparam T Type to construct
  /// @tparam Args Constructor argument types
  /// @param args Arguments to forward to T's constructor
  /// @return Reference to the newly constructed value
  template <typename T, typename... Args>
    requires(std::is_copy_assignable_v<std::decay_t<T>> &&
             std::is_copy_constructible_v<std::decay_t<T>> &&
             (sizeof(std::decay_t<T>) > SIZE ||
              (std::is_move_assignable_v<std::decay_t<T>> &&
               std::is_move_constructible_v<std::decay_t<T>>)))
  T& emplace(Args&&... args) {
    using U = std::decay_t<T>;
    destroy();
    m_handler = makeHandler<U>();
    return *constructValue<U>(std::forward<Args>(args)...);
  }

  /// Get reference to stored value of specified type
  /// @tparam T Type to retrieve (must be exact type, no const/ref)
  /// @return Reference to the stored value
  /// @throws std::bad_any_cast if stored type doesn't match T
  template <typename T>
  T& as() {
    static_assert(std::is_same_v<T, std::decay_t<T>>,
                  "Please pass the raw type, no const or ref");
    if (m_handler == nullptr || m_handler->typeHash != typeHash<T>()) {
      throw std::bad_any_cast{};
    }

    _ACTS_ANY_VERBOSE("Get as "
                      << (m_handler->heapAllocated ? "heap" : "local"));

    return *reinterpret_cast<T*>(dataPtr());
  }

  /// Get const reference to stored value of specified type
  /// @tparam T Type to retrieve (must be exact type, no const/ref)
  /// @return Const reference to the stored value
  /// @throws std::bad_any_cast if stored type doesn't match T
  template <typename T>
  const T& as() const {
    static_assert(std::is_same_v<T, std::decay_t<T>>,
                  "Please pass the raw type, no const or ref");
    if (m_handler == nullptr || m_handler->typeHash != typeHash<T>()) {
      throw std::bad_any_cast{};
    }

    _ACTS_ANY_VERBOSE("Get as " << (m_handler->heap ? "heap" : "local"));

    return *reinterpret_cast<const T*>(dataPtr());
  }

  ~AnyBase() { destroy(); }

  /// Copy constructor
  /// @param other The AnyBase to copy from
  AnyBase(const AnyBase& other) _ACTS_ANY_NOEXCEPT {
    if (m_handler == nullptr && other.m_handler == nullptr) {
      // both are empty, noop
      return;
    }

    _ACTS_ANY_VERBOSE("Copy construct (this="
                      << this << ") at: " << static_cast<void*>(m_data.data()));

    m_handler = other.m_handler;
    copyConstruct(other);
  }

  /// Copy assignment operator
  /// @param other The AnyBase to copy from
  /// @return Reference to this object
  AnyBase& operator=(const AnyBase& other) _ACTS_ANY_NOEXCEPT {
    _ACTS_ANY_VERBOSE("Copy assign (this="
                      << this << ") at: " << static_cast<void*>(m_data.data()));

    if (m_handler == nullptr && other.m_handler == nullptr) {
      // both are empty, noop
      return *this;
    }

    if (m_handler == other.m_handler) {
      // same type, but checked before they're not both nullptr
      copy(std::move(other));
    } else {
      if (m_handler != nullptr) {
        // this object is not empty, but have different types => destroy
        destroy();
      }
      assert(m_handler == nullptr);
      m_handler = other.m_handler;
      copyConstruct(std::move(other));
    }
    return *this;
  }

  /// Move constructor
  /// @param other The AnyBase to move from
  AnyBase(AnyBase&& other) _ACTS_ANY_NOEXCEPT {
    _ACTS_ANY_VERBOSE("Move construct (this="
                      << this << ") at: " << static_cast<void*>(m_data.data()));
    if (m_handler == nullptr && other.m_handler == nullptr) {
      // both are empty, noop
      return;
    }

    m_handler = other.m_handler;
    moveConstruct(std::move(other));
  }

  /// Move assignment operator
  /// @param other The AnyBase to move from
  /// @return Reference to this object
  AnyBase& operator=(AnyBase&& other) _ACTS_ANY_NOEXCEPT {
    _ACTS_ANY_VERBOSE("Move assign (this="
                      << this << ") at: " << static_cast<void*>(m_data.data()));
    if (m_handler == nullptr && other.m_handler == nullptr) {
      // both are empty, noop
      return *this;
    }

    if (m_handler == other.m_handler) {
      // same type, but checked before they're not both nullptr
      move(std::move(other));
    } else {
      if (m_handler != nullptr) {
        // this object is not empty, but have different types => destroy
        destroy();
      }
      assert(m_handler == nullptr);
      m_handler = other.m_handler;
      moveConstruct(std::move(other));
    }

    return *this;
  }

  /// Check if the AnyBase contains a value
  /// @return True if a value is stored, false if empty
  explicit operator bool() const { return m_handler != nullptr; }

 private:
  void* dataPtr() {
    if (m_handler->heapAllocated) {
      return *reinterpret_cast<void**>(m_data.data());
    } else {
      return reinterpret_cast<void*>(m_data.data());
    }
  }

  void setDataPtr(void* ptr) { *reinterpret_cast<void**>(m_data.data()) = ptr; }

  const void* dataPtr() const {
    if (m_handler->heapAllocated) {
      return *reinterpret_cast<void* const*>(m_data.data());
    } else {
      return reinterpret_cast<const void*>(m_data.data());
    }
  }

  template <typename T>
  static std::uint64_t typeHash() {
    const static std::uint64_t value = detail::fnv1a_64(typeid(T).name());
    return value;
  }

  struct Handler {
    void (*destroy)(void* ptr) = nullptr;
    void (*moveConstruct)(void* from, void* to) = nullptr;
    void (*move)(void* from, void* to) = nullptr;
    void* (*copyConstruct)(const void* from, void* to) = nullptr;
    void (*copy)(const void* from, void* to) = nullptr;
    bool heapAllocated{false};
    std::uint64_t typeHash{0};
  };

  template <typename T>
  static const Handler* makeHandler() {
    static_assert(!std::is_same_v<T, AnyBase<SIZE>>, "Cannot wrap any in any");
    static const Handler static_handler = []() {
      Handler h;
      h.heapAllocated = heapAllocated<T>();
      if constexpr (!std::is_trivially_destructible_v<T> ||
                    heapAllocated<T>()) {
        h.destroy = &destroyImpl<T>;
      }
      if constexpr (!heapAllocated<T>() &&
                    !std::is_trivially_move_constructible_v<T>) {
        h.moveConstruct = &moveConstructImpl<T>;
      }
      if constexpr (!heapAllocated<T>() &&
                    !std::is_trivially_move_assignable_v<T>) {
        h.move = &moveImpl<T>;
      }
      if constexpr (!std::is_trivially_copy_constructible_v<T> ||
                    heapAllocated<T>()) {
        h.copyConstruct = &copyConstructImpl<T>;
      }
      if constexpr (!std::is_trivially_copy_assignable_v<T> ||
                    heapAllocated<T>()) {
        h.copy = &copyImpl<T>;
      }

      h.typeHash = typeHash<T>();

      _ACTS_ANY_DEBUG("Type: " << typeid(T).name());
      _ACTS_ANY_DEBUG(" -> destroy: " << h.destroy);
      _ACTS_ANY_DEBUG(" -> moveConstruct: " << h.moveConstruct);
      _ACTS_ANY_DEBUG(" -> move: " << h.move);
      _ACTS_ANY_DEBUG(" -> copyConstruct: " << h.copyConstruct);
      _ACTS_ANY_DEBUG(" -> copy: " << h.copy);
      _ACTS_ANY_DEBUG(
          " -> heapAllocated: " << (h.heapAllocated ? "yes" : "no"));

      return h;
    }();
    return &static_handler;
  }

  template <typename T>
  static constexpr bool heapAllocated() {
    return sizeof(T) > SIZE;
  }

  template <typename T, typename... Args>
  T* constructValue(Args&&... args) {
    if constexpr (!heapAllocated<T>()) {
      // construct into local buffer
      auto* ptr = new (m_data.data()) T(std::forward<Args>(args)...);
      _ACTS_ANY_VERBOSE("Construct local (this="
                        << this
                        << ") at: " << static_cast<void*>(m_data.data()));
      return ptr;
    } else {
      // too large, heap allocate
      T* heap = new T(std::forward<Args>(args)...);
      _ACTS_ANY_DEBUG("Allocate type: " << typeid(T).name() << " at " << heap);
      _ACTS_ANY_TRACK_ALLOCATION(T, heap);
      setDataPtr(heap);
      return heap;
    }
  }

  void destroy() {
    _ACTS_ANY_VERBOSE("Destructor this=" << this << " handler: " << m_handler);
    if (m_handler != nullptr && m_handler->destroy != nullptr) {
      _ACTS_ANY_VERBOSE("Non-trivial destruction");
      m_handler->destroy(dataPtr());
    }
    m_handler = nullptr;
  }

  void moveConstruct(AnyBase&& fromAny) {
    if (m_handler == nullptr) {
      return;
    }

    void* to = dataPtr();
    void* from = fromAny.dataPtr();
    if (m_handler->heapAllocated) {
      // stored on heap: just copy the pointer
      setDataPtr(fromAny.dataPtr());
      // do not delete in moved-from any
      fromAny.m_handler = nullptr;
      return;
    }

    if (m_handler->moveConstruct == nullptr) {
      _ACTS_ANY_VERBOSE("Trivially move construct");
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

    void* to = dataPtr();
    void* from = fromAny.dataPtr();
    if (m_handler->heapAllocated) {
      // stored on heap: just copy the pointer
      // need to delete existing pointer
      m_handler->destroy(dataPtr());
      setDataPtr(fromAny.dataPtr());
      // do not delete in moved-from any
      fromAny.m_handler = nullptr;
      return;
    }

    if (m_handler->move == nullptr) {
      _ACTS_ANY_VERBOSE("Trivially move");
      // trivially movable
      m_data = std::move(fromAny.m_data);
    } else {
      m_handler->move(from, to);
    }
  }

  void copyConstruct(const AnyBase& fromAny) {
    if (m_handler == nullptr) {
      return;
    }

    void* to = dataPtr();
    const void* from = fromAny.dataPtr();

    if (m_handler->copyConstruct == nullptr) {
      _ACTS_ANY_VERBOSE("Trivially copy construct");
      // trivially copy constructible
      m_data = fromAny.m_data;
    } else {
      void* copyAt = m_handler->copyConstruct(from, to);
      if (to == nullptr) {
        assert(copyAt != nullptr);
        // copy allocated, store pointer
        setDataPtr(copyAt);
      }
    }
  }

  void copy(const AnyBase& fromAny) {
    if (m_handler == nullptr) {
      return;
    }

    void* to = dataPtr();
    const void* from = fromAny.dataPtr();

    if (m_handler->copy == nullptr) {
      _ACTS_ANY_VERBOSE("Trivially copy");
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
    if constexpr (!heapAllocated<T>()) {
      // stored in place: just call the destructor
      _ACTS_ANY_VERBOSE("Destroy local at: " << ptr);
      obj->~T();
    } else {
      // stored on heap: delete
      _ACTS_ANY_DEBUG("Delete type: " << typeid(T).name()
                                      << " heap at: " << obj);
      _ACTS_ANY_TRACK_DEALLOCATION(T, obj);
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
    const T* _from = static_cast<const T*>(from);
    if (to == nullptr) {
      assert(heapAllocated<T>() && "Received nullptr in local buffer case");
      to = new T(*_from);
      _ACTS_ANY_TRACK_ALLOCATION(T, to);

    } else {
      assert(!heapAllocated<T>() && "Received non-nullptr in heap case");
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

  static constexpr std::size_t kMaxAlignment =
      std::max(alignof(std::max_align_t),
#if defined(__AVX512F__)
               std::size_t{64}
#elif defined(__AVX__)
               std::size_t{32}
#elif defined(__SSE__)
               std::size_t{16}
#else
               std::size_t{0}
  // Neutral element
  // for maximum
#endif
      );

  alignas(kMaxAlignment) std::array<std::byte, SIZE> m_data{};
  const Handler* m_handler{nullptr};
};

/// @brief A type-safe container for single values of any type
/// @details This is a custom implementation similar to `std::any` but optimized for small types
///          that can fit into a pointer-sized buffer. Values larger than a
///          pointer are stored on the heap.
using Any = AnyBase<sizeof(void*)>;

/// @}

#undef _ACTS_ANY_VERBOSE
#undef _ACTS_ANY_VERBOSE_BUFFER
#undef _ACTS_ANY_ENABLE_VERBOSE

}  // namespace Acts

// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Utilities/TypeTraits.hpp"

#include <cassert>
#include <functional>
#include <memory>
#include <type_traits>

namespace Acts {

/// Ownership enum for @c Delegate
enum class DelegateType { Owning, NonOwning };

// Specialization needed for defaulting ownership and for R(Args...) syntax
template <typename, typename H = void, DelegateType = DelegateType::NonOwning>
class Delegate;

/// Delegate type that allows type erasure of a callable without allocation
/// and with a single level of indirection. This type can support:
/// - a free function pointer
/// - a pointer to a member function alongside an instance pointer
/// @note @c Delegate by default does not assume ownership of the instance.
///          In that case You need to ensure that the lifetime of the callable
///          instance is longer than that of the @c Delegate. If you set @c O
///          to @c DelegateType::Owning, it will assume ownership.
/// @note Currently @c Delegate only supports callables that are ``const``
/// @tparam R Return type of the function signature
/// @tparam H Holder type that is used to store an instance
/// @tparam O Ownership type of the delegate: Owning or NonOwning
/// @tparam Args Types of the arguments of the function signatures
///
template <typename R, typename H, DelegateType O, typename... Args>
class Delegate<R(Args...), H, O> {
  static constexpr DelegateType kOwnership = O;

  /// Alias of the return type
  using return_type = R;
  using holder_type = H;
  /// Alias to the function pointer type this class will store
  using function_type = return_type (*)(const holder_type *, Args...);

  using function_ptr_type = return_type (*)(Args...);

  using deleter_type = void (*)(const holder_type *);

  template <typename T, typename C>
  using isSignatureCompatible =
      decltype(std::declval<T &>() = std::declval<C>());

  using OwningDelegate =
      Delegate<R(Args...), holder_type, DelegateType::Owning>;
  using NonOwningDelegate =
      Delegate<R(Args...), holder_type, DelegateType::NonOwning>;
  template <typename T>
  using isNoFunPtr = std::enable_if_t<
      not std::is_convertible_v<std::decay_t<T>, function_type> and
      not std::is_same_v<std::decay_t<T>, OwningDelegate> and
      not std::is_same_v<std::decay_t<T>, NonOwningDelegate>>;

 public:
  Delegate() = default;

  Delegate(Delegate &&) = default;
  Delegate &operator=(Delegate &&) = default;
  Delegate(const Delegate &) = default;
  Delegate &operator=(const Delegate &) = default;

  /// Constructor with an explicit runtime callable
  /// @param callable The runtime value of the callable
  /// @note The function signature requires the first argument of the callable is `const void*`.
  ///       i.e. if the signature of the delegate is `void(int)`, the
  ///       callable's signature has to be `void(const void*, int)`.
  Delegate(function_type callable) { connect(callable); }

  /// Constructor with a possibly stateful function object.
  /// @tparam Callable Type of the callable
  /// @param callable The callable (function object or lambda)
  /// @note @c Delegate does not assume owner ship over @p callable. You need to ensure
  ///       it's lifetime is longer than that of @c Delegate.
  template <typename Callable, typename = isNoFunPtr<Callable>>
  Delegate(Callable &callable) {
    connect(callable);
  }

  /// Constructor from rvalue reference is deleted, should catch construction
  /// with temporary objects and thus invalid pointers
  template <typename Callable, typename = isNoFunPtr<Callable>>
  Delegate(Callable &&) = delete;

  /// Assignment operator with an explicit runtime callable
  /// @param callable The runtime value of the callable
  /// @note The function signature requires the first argument of the callable is `const void*`.
  ///       i.e. if the signature of the delegate is `void(int)`, the
  ///       callable's signature has to be `void(const void*, int)`.
  void operator=(function_type callable) { connect(callable); }

  /// Assignment operator with possibly stateful function object.
  /// @tparam Callable Type of the callable
  /// @param callable The callable (function object or lambda)
  /// @note @c Delegate does not assume owner ship over @p callable. You need to ensure
  ///       it's lifetime is longer than that of @c Delegate.
  template <typename Callable, typename = isNoFunPtr<Callable>>
  void operator=(Callable &callable) {
    connect(callable);
  }

  /// Assignment operator from rvalue reference is deleted, should catch
  /// assignment from temporary objects and thus invalid pointers
  template <typename Callable, typename = isNoFunPtr<Callable>>
  void operator=(Callable &&) = delete;

  /// Connect a free function pointer.
  /// @note The function pointer must be ``constexpr`` for @c Delegate to accept it
  /// @tparam Callable The compile-time free function pointer
  template <auto Callable>
  void connect() {
    m_payload.payload = nullptr;

    static_assert(
        Concepts::is_detected<isSignatureCompatible, function_ptr_type,
                              decltype(Callable)>::value,
        "Callable given does not correspond exactly to required call "
        "signature");

    m_function = [](const holder_type * /*payload*/,
                    Args... args) -> return_type {
      return std::invoke(Callable, std::forward<Args>(args)...);
    };
  }

  /// Assignment operator with possibly stateful function object.
  /// @tparam Callable Type of the callable
  /// @param callable The callable (function object or lambda)
  /// @note @c Delegate does not assume owner ship over @p callable. You need to ensure
  ///       it's lifetime is longer than that of @c Delegate.
  template <typename Callable, typename = isNoFunPtr<Callable>>
  void connect(Callable &callable) {
    connect<&Callable::operator(), Callable>(&callable);
  }

  /// Connection with rvalue reference is deleted, should catch assignment
  /// from temporary objects and thus invalid pointers
  template <typename Callable, typename = isNoFunPtr<Callable>>
  void connect(Callable &&) = delete;

  /// Connect anything that is assignable to the function pointer
  /// @param callable The runtime value of the callable
  /// @note The function signature requires the first argument of the callable is `const void*`.
  ///       i.e. if the signature of the delegate is `void(int)`, the
  ///       callable's signature has to be `void(const void*, int)`.
  void connect(function_type callable) {
    if constexpr (kOwnership == DelegateType::NonOwning) {
      m_payload.payload = nullptr;
    }
    m_function = callable;
  }

  /// Connect a member function to be called on an instance
  /// @tparam Callable The compile-time member function pointer
  /// @tparam Type The type of the instance the member function should be called on
  /// @param instance The instance on which the member function pointer should be called on
  /// @note @c Delegate does not assume owner ship over @p instance. You need to ensure
  ///       it's lifetime is longer than that of @c Delegate.
  template <auto Callable, typename Type, DelegateType T = kOwnership,
            typename = std::enable_if_t<T == DelegateType::NonOwning>>
  void connect(const Type *instance) {
    using member_ptr_type = return_type (Type::*)(Args...) const;

    static_assert(Concepts::is_detected<isSignatureCompatible, member_ptr_type,
                                        decltype(Callable)>::value,
                  "Callable given does not correspond exactly to required call "
                  "signature");

    m_payload.payload = instance;

    m_function = [](const holder_type *payload, Args... args) -> return_type {
      assert(payload != nullptr && "Payload is required, but not set");
      const auto *concretePayload = static_cast<const Type *>(payload);
      return std::invoke(Callable, concretePayload,
                         std::forward<Args>(args)...);
    };
  }

  /// Connect a member function to be called on an instance
  /// @tparam Callable The compile-time member function pointer
  /// @tparam Type The type of the instance the member function should be called on
  /// @param instance The instance on which the member function pointer should be called on
  /// @note @c Delegate assumes owner ship over @p instance.
  template <auto Callable, typename Type>
  void connect(std::unique_ptr<const Type> instance) {
    using member_ptr_type = return_type (Type::*)(Args...) const;
    static_assert(Concepts::is_detected<isSignatureCompatible, member_ptr_type,
                                        decltype(Callable)>::value,
                  "Callable given does not correspond exactly to required call "
                  "signature");

    if constexpr (kOwnership == DelegateType::Owning) {
      m_payload.payload = std::unique_ptr<const holder_type, deleter_type>(
          instance.release(), [](const holder_type *payload) {
            const auto *concretePayload = static_cast<const Type *>(payload);
            delete concretePayload;
          });
    } else {
      m_payload.payload = instance.release();
    }

    m_function = [](const holder_type *payload, Args... args) -> return_type {
      assert(payload != nullptr && "Payload is required, but not set");
      const auto *concretePayload = static_cast<const Type *>(payload);
      return std::invoke(Callable, concretePayload,
                         std::forward<Args>(args)...);
    };
  }

  /// The call operator that exposes the functionality of the @c Delegate type.
  /// @param args The arguments to call the contained function with
  /// @return Return value of the contained function
  return_type operator()(Args... args) const {
    assert(connected() && "Delegate is not connected");
    return std::invoke(m_function, m_payload.ptr(),
                       std::forward<Args>(args)...);
  }

  /// Return whether this delegate is currently connected
  /// @return True if this delegate is connected
  bool connected() const { return m_function != nullptr; }

  /// Return whether this delegate is currently connected
  /// @return True if this delegate is connected
  operator bool() const { return connected(); }

  /// Disconnect this delegate, meaning it cannot be called anymore
  void disconnect() {
    m_payload.clear();
    m_function = nullptr;
  }

  template <typename holder_t = holder_type,
            typename = std::enable_if_t<!std::is_same_v<holder_t, void>>>
  const holder_type *instance() const {
    return m_payload.ptr();
  }

 private:
  // Deleter that does not do anything
  static void noopDeleter(const holder_type * /*unused*/) {}

  /// @cond

  // Payload object without a deleter
  struct NonOwningPayload {
    void clear() { payload = nullptr; }

    const holder_type *ptr() const { return payload; }

    const holder_type *payload{nullptr};
  };

  // Payload object with a deleter
  struct OwningPayload {
    void clear() { payload.reset(); }

    const holder_type *ptr() const { return payload.get(); }

    std::unique_ptr<const holder_type, deleter_type> payload{nullptr,
                                                             &noopDeleter};
  };

  /// Stores the instance pointer and maybe a deleter
  std::conditional_t<kOwnership == DelegateType::NonOwning, NonOwningPayload,
                     OwningPayload>
      m_payload;

  /// @endcond

  /// Stores the function pointer wrapping the compile time function pointer given in @c connect().
  function_type m_function{nullptr};
};

template <typename, typename H = void>
class OwningDelegate;

/// Alias for an owning delegate
template <typename R, typename H, typename... Args>
class OwningDelegate<R(Args...), H>
    : public Delegate<R(Args...), H, DelegateType::Owning> {};

}  // namespace Acts

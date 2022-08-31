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

namespace Acts {

template <typename>
class Delegate;

/// Delegate type that allows type erasure of a callable without allocation and
/// with a single level of indirection.
/// This type can support:
/// - a free function pointer
/// - a pointer to a member function alongside an instance pointer
/// @note @c Delegate does not assume ownership of the instance.
///          You need to ensure that the lifetime of the callable
///          instance is longer than that of the @c Delegate.
/// @note Currently @c Delegate only supports callables that are ``const``
/// @tparam R Return type of the function signature
/// @tparam Args Types of the arguments of the function signatures
///
template <typename R, typename... Args>
class Delegate<R(Args...)> {
  /// Alias of the return type
  using return_type = R;
  /// Alias to the function pointer type this class will store
  using function_type = return_type (*)(const void *, Args...);

  using function_ptr_type = return_type (*)(Args...);

  template <typename T, typename C>
  using isSignatureCompatible =
      decltype(std::declval<T &>() = std::declval<C>());

  template <typename T>
  using isNoFunPtr = std::enable_if_t<
      not std::is_convertible_v<std::decay_t<T>, function_type>>;

 public:
  Delegate() = default;

  /// Constructor with an explicit runtime callable
  /// @param callable The runtime value of the callable
  /// @note The function signature requires the first argument of the callable is `const void*`.
  ///       i.e. if the signature of the delegate is `void(int)`, the callable's
  ///       signature has to be `void(const void*, int)`.
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
  ///       i.e. if the signature of the delegate is `void(int)`, the callable's
  ///       signature has to be `void(const void*, int)`.
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
    m_payload = nullptr;

    static_assert(
        Concepts::is_detected<isSignatureCompatible, function_ptr_type,
                              decltype(Callable)>::value,
        "Callable given does not correspond exactly to required call "
        "signature");

    m_function = [](const void * /*payload*/, Args... args) -> return_type {
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

  /// Connection with rvalue reference is deleted, should catch assignment from
  /// temporary objects and thus invalid pointers
  template <typename Callable, typename = isNoFunPtr<Callable>>
  void connect(Callable &&) = delete;

  /// Connect anything that is assignable to the function pointer
  /// @param callable The runtime value of the callable
  /// @note The function signature requires the first argument of the callable is `const void*`.
  ///       i.e. if the signature of the delegate is `void(int)`, the callable's
  ///       signature has to be `void(const void*, int)`.
  void connect(function_type callable) {
    m_payload = nullptr;
    m_function = callable;
  }

  /// Connect a member function to be called on an instance
  /// @tparam Callable The compile-time member function pointer
  /// @tparam Type The type of the instance the member function should be called on
  /// @param instance The instance on which the member function pointer should be called on
  /// @note @c Delegate does not assume owner ship over @p instance. You need to ensure
  ///       it's lifetime is longer than that of @c Delegate.
  template <auto Callable, typename Type>
  void connect(const Type *instance) {
    using member_ptr_type = return_type (Type::*)(Args...) const;

    static_assert(Concepts::is_detected<isSignatureCompatible, member_ptr_type,
                                        decltype(Callable)>::value,
                  "Callable given does not correspond exactly to required call "
                  "signature");

    m_payload = instance;
    m_function = [](const void *payload, Args... args) -> return_type {
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
    return std::invoke(m_function, m_payload, std::forward<Args>(args)...);
  }

  /// Return whether this delegate is currently connected
  /// @return True if this delegate is connected
  bool connected() const { return m_function != nullptr; }

  /// Return whether this delegate is currently connected
  /// @return True if this delegate is connected
  operator bool() const { return connected(); }

  /// Disconnect this delegate, meaning it cannot be called anymore
  void disconnect() {
    m_payload = nullptr;
    m_function = nullptr;
  }

 private:
  /// Stores the instance pointer
  const void *m_payload{nullptr};
  /// Stores the function pointer wrapping the compile time function pointer given in @c connect().
  function_type m_function{nullptr};
};
}  // namespace Acts

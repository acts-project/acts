// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include <boost/core/demangle.hpp>
#include <format>
#include <functional>
#include <typeindex>
#include <typeinfo>
#include <unordered_map>

namespace Acts {

/// Template class for type-based function dispatch
/// 
/// This class allows registering function pointers associated with specific
/// derived types of a base class. When invoked with a base class reference,
/// it will look up and call the appropriate registered function.
/// 
/// @tparam BaseType The base class type that will be the first parameter
/// @tparam Signature The function signature (e.g., ReturnType(Args...))
template <typename BaseType, typename Signature>
class TypeDispatcher;

template <typename BaseType, typename ReturnType, typename... Args>
class TypeDispatcher<BaseType, ReturnType(Args...)> {
 public:
  using function_signature = ReturnType(const BaseType&, Args...);
  using function_type = std::function<function_signature>;

  /// Register a free function for a specific derived type
  /// @tparam DerivedType The derived type to associate the function with
  /// @param func The function to register
  template <typename DerivedType>
    requires std::is_base_of_v<BaseType, DerivedType>
  void registerFunction(ReturnType (*func)(const DerivedType&, Args...)) {
    registerFunctionImpl<DerivedType>([func](const DerivedType& derived, Args... args) -> ReturnType {
      return func(derived, std::forward<Args>(args)...);
    });
  }

  /// Register a function object or lambda for a specific derived type (by rvalue)
  /// @tparam DerivedType The derived type to associate the function with
  /// @tparam FuncType The function object type
  /// @param func The function object to register (must be rvalue to indicate ownership transfer)
  template <typename DerivedType, typename FuncType>
    requires std::is_base_of_v<BaseType, DerivedType>
  void registerFunction(FuncType&& func) {
    registerFunctionImpl<DerivedType>(std::move(func));
  }

  /// Register a std::function for a specific derived type (by rvalue only)
  /// @tparam DerivedType The derived type to associate the function with
  /// @param func The std::function to register (must be rvalue to indicate ownership transfer)
  template <typename DerivedType>
    requires std::is_base_of_v<BaseType, DerivedType>
  void registerFunction(std::function<ReturnType(const DerivedType&, Args...)>&& func) {
    registerFunctionImpl<DerivedType>(std::move(func));
  }

  /// Explicitly deleted overload for lvalue std::function to enforce rvalue-only registration.
  /// Without this, lvalue std::functions would match the universal reference overload above,
  /// which would copy the function instead of moving it, defeating our ownership semantics.
  /// This provides a clear compilation error with a helpful message instead.
  template <typename DerivedType>
  void registerFunction(const std::function<ReturnType(const DerivedType&, Args...)>& func) = delete;

  /// Call the registered function for the given object's type
  /// @param obj The object to dispatch on
  /// @param args Additional arguments to pass to the function
  /// @return The return value from the registered function
  ReturnType operator()(const BaseType& obj, Args... args) const {
    std::type_index typeIdx(typeid(obj));
    
    auto it = m_functions.find(typeIdx);
    if (it == m_functions.end()) {
      throw std::runtime_error(std::format("No function registered for type: {}", 
                                           boost::core::demangle(typeIdx.name())));
    }
    
    return it->second(obj, std::forward<Args>(args)...);
  }

  /// Check if a function is registered for the given object's type
  /// @param obj The object to check
  /// @return true if a function is registered, false otherwise
  bool hasFunction(const BaseType& obj) const {
    std::type_index typeIdx(typeid(obj));
    return m_functions.find(typeIdx) != m_functions.end();
  }

  /// Check if a function is registered for the given type
  /// @tparam DerivedType The type to check
  /// @return true if a function is registered, false otherwise
  template <typename DerivedType>
    requires std::is_base_of_v<BaseType, DerivedType>
  bool hasFunction() const {
    std::type_index typeIdx(typeid(DerivedType));
    return m_functions.find(typeIdx) != m_functions.end();
  }

  /// Clear all registered functions
  void clear() {
    m_functions.clear();
  }

  /// Get the number of registered functions
  std::size_t size() const {
    return m_functions.size();
  }

 private:
  /// Internal helper method that centralizes function registration logic
  /// @tparam DerivedType The derived type to associate the function with
  /// @tparam FuncType The function type (deduced)
  /// @param func The function to register
  template <typename DerivedType, typename FuncType>
  void registerFunctionImpl(FuncType&& func) {
    std::type_index typeIdx(typeid(DerivedType));
    
    // Wrap the function in a lambda that performs the dynamic cast
    m_functions[typeIdx] = [func = std::forward<FuncType>(func)](
        const BaseType& base, Args... args) -> ReturnType {
      const auto* derived = dynamic_cast<const DerivedType*>(&base);
      if (derived == nullptr) {
        throw std::bad_cast();
      }
      return func(*derived, std::forward<Args>(args)...);
    };
  }

  std::unordered_map<std::type_index, function_type> m_functions;
};

}  // namespace Acts
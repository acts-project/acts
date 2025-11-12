// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include <format>
#include <functional>
#include <type_traits>
#include <typeindex>
#include <typeinfo>
#include <unordered_map>
#include <vector>

#include <boost/core/demangle.hpp>

namespace Acts {

/// Template class for type-based function dispatch
///
/// This class allows registering function pointers associated with specific
/// derived types of a base class. When invoked with a base class reference,
/// it will look up and call the appropriate registered function.
///
/// @tparam base_t The base class type that will be the first parameter
/// @tparam signature_t The function signature (e.g., return_t(Args...))
template <typename base_t, typename signature_t>
class TypeDispatcher;

template <typename base_t, typename return_t, typename... args_t>
class TypeDispatcher<base_t, return_t(args_t...)> {
 public:
  // Type aliases for frequently used template parameters
  using base_type = base_t;
  using return_type = return_t;
  using self_type = TypeDispatcher<base_type, return_type(args_t...)>;

  using function_signature = return_t(const base_t&, args_t...);
  using function_pointer_type = return_t (*)(args_t...);
  using function_type = std::function<function_signature>;

  /// Default constructor
  TypeDispatcher() = default;

  /// Constructor that registers multiple function pointers with auto-detected
  /// types
  template <typename... derived_types_t>
  explicit TypeDispatcher(return_t (*... funcs)(const derived_types_t&,
                                                args_t...))
    requires(std::is_base_of_v<base_t, derived_types_t> && ...)
  {
    (registerFunction(funcs),
     ...);  // Fold expression to register each function
  }

  /// Register a free function with explicit derived type
  /// @tparam derived_t The derived type to associate the function with
  /// @param func The function pointer
  template <typename derived_t>
    requires std::is_base_of_v<base_t, derived_t>
  self_type& registerFunction(return_t (*func)(const derived_t&, args_t...)) {
    std::type_index typeIdx(typeid(derived_t));

    // Check if this exact type is already registered
    if (m_functions.find(typeIdx) != m_functions.end()) {
      throw std::runtime_error(
          std::format("Function already registered for type: {}",
                      boost::core::demangle(typeIdx.name())));
    }

    // Try to detect conflicts with existing registrations if the type is
    // default constructible
    if constexpr (std::is_default_constructible_v<derived_t>) {
      derived_t tempObj{};
      for (const auto& [existingTypeIdx, checker] : m_castCheckers) {
        if (checker(tempObj)) {
          throw std::runtime_error(
              std::format("Registration conflict: type {} would be handled by "
                          "existing function registered for type: {}",
                          boost::core::demangle(typeIdx.name()),
                          boost::core::demangle(existingTypeIdx.name())));
        }
      }
    }

    // Store a cast checker that tests if dynamic_cast<derived_t*> will work
    m_castCheckers[typeIdx] = [](const base_t& obj) -> bool {
      const auto* derived = dynamic_cast<const derived_t*>(&obj);
      return derived != nullptr;
    };

    // Wrap the function in a lambda that performs the dynamic cast
    m_functions[typeIdx] = [func]<typename... Ts>(const base_t& base,
                                                  Ts&&... args) -> return_t {
      const auto* derived = dynamic_cast<const derived_t*>(&base);
      if (derived == nullptr) {
        throw std::bad_cast();
      }
      return func(*derived, std::forward<Ts>(args)...);
    };

    return *this;
  }

  /// Call the registered function for the given object's type
  /// @param obj The object to dispatch on
  /// @param args Additional arguments to pass to the function
  /// @return The return value from the registered function
  template <typename... func_args_t>
  return_t operator()(const base_t& obj, func_args_t&&... args) const
    requires std::invocable<function_pointer_type, func_args_t...>
  {
    std::vector<std::type_index> compatibleTypes;

    // Find all registered functions that can handle this object type
    for (const auto& [registeredTypeIdx, checker] : m_castCheckers) {
      if (checker(obj)) {
        compatibleTypes.push_back(registeredTypeIdx);
      }
    }

    if (compatibleTypes.empty()) {
      throw std::runtime_error(
          std::format("No function registered for type: {}",
                      boost::core::demangle(typeid(obj).name())));
    }

    if (compatibleTypes.size() > 1) {
      std::string typeNames;
      for (std::size_t i = 0; i < compatibleTypes.size(); ++i) {
        if (i > 0) {
          typeNames += ", ";
        }
        typeNames += boost::core::demangle(compatibleTypes[i].name());
      }
      throw std::runtime_error(
          std::format("Ambiguous dispatch for type {}: multiple functions can "
                      "handle it: {}",
                      boost::core::demangle(typeid(obj).name()), typeNames));
    }

    // Exactly one compatible function found
    auto funcIt = m_functions.find(compatibleTypes[0]);
    if (funcIt != m_functions.end()) {
      return funcIt->second(obj, std::forward<func_args_t>(args)...);
    }

    // This should never happen if our data structures are consistent
    throw std::runtime_error(
        "Internal error: function not found for compatible type");
  }

  /// Check if a function is registered for the given object's type
  /// @param obj The object to check
  /// @return true if a function is registered, false otherwise
  bool hasFunction(const base_t& obj) const {
    // Find all registered functions that can handle this object type
    for (const auto& [registeredTypeIdx, checker] : m_castCheckers) {
      if (checker(obj)) {
        return true;
      }
    }
    return false;
  }

  /// Check if a function is registered for the given type
  /// @tparam derived_t The type to check
  /// @return true if a function is registered, false otherwise
  template <typename derived_t>
    requires std::is_base_of_v<base_t, derived_t>
  bool hasFunction() const {
    std::type_index typeIdx(typeid(derived_t));
    return m_functions.find(typeIdx) != m_functions.end();
  }

  /// Clear all registered functions
  void clear() {
    m_functions.clear();
    m_castCheckers.clear();
  }

  /// Get the number of registered functions
  std::size_t size() const { return m_functions.size(); }

 private:
  using cast_checker_type = std::function<bool(const base_t&)>;

  std::unordered_map<std::type_index, function_type> m_functions;
  std::unordered_map<std::type_index, cast_checker_type> m_castCheckers;
};

}  // namespace Acts

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
#include <optional>
#include <stdexcept>
#include <type_traits>
#include <typeindex>
#include <typeinfo>
#include <utility>
#include <vector>

#include <boost/core/demangle.hpp>

namespace Acts {

/// Template class for symmetric dispatch on a *pair* of polymorphic types.
///
/// This is a sibling of @c TypeDispatcher. Where @c TypeDispatcher selects a
/// function based on the concrete type of a single base-class reference, this
/// class selects a function based on the concrete types of *two* base-class
/// references.
///
/// Registered functions are treated as **symmetric**: a function registered
/// with parameter types @c (A, B) will be selected both when called as
/// @c (a, b) and as @c (b, a). In either case the registered function is
/// invoked with its arguments in the order it declared them, i.e. the @c A
/// argument first and the @c B argument second. This lets callers register a
/// single handler per unordered type pair and rely on the dispatcher to
/// normalise the argument order.
///
/// @note Matching is based on @c dynamic_cast, so a handler parameter typed as
/// a base (the root base or an intermediate class) matches *any* object in that
/// slot: a @c dynamic_cast to a base always succeeds. A base-typed slot is
/// therefore a **wildcard**, and @c (Base, B) behaves as "a @c B paired with
/// anything".
///
/// @note There is **no specificity ranking** — this dispatcher does not mirror
/// C++ overload resolution and there is no "most derived wins" tie-break. All
/// matching registrations are considered equal; if more than one distinct
/// registration matches a call it is reported as an ambiguity rather than being
/// resolved in favour of the more derived handler. Standard RTTI exposes no
/// base-class information, so the required subtype order cannot be recovered
/// from the registered types in the general case. Consequently, callers should
/// **register non-overlapping type pairs** and avoid mixing a wildcard
/// (base-typed) slot with more specific handlers for the same subtree.
///
/// Overlaps are best-effort rejected at registration time, but only when the
/// newly registered derived types are default constructible (a sample of each
/// is materialised to probe existing registrations). Overlaps that cannot be
/// detected this way — e.g. involving an abstract base or a non-default-
/// constructible type, or where the wildcard handler is registered *last* —
/// surface only at call time as an ambiguity, and only for the specific pairs
/// that actually trigger it.
///
/// @tparam base_t The common base class of both dispatch arguments
/// @tparam signature_t The trailing function signature (e.g. return_t(Args...))
template <typename base_t, typename signature_t>
class PairTypeDispatcher;

/// Pair type dispatcher specialization for function signature
template <typename base_t, typename return_t, typename... args_t>
class PairTypeDispatcher<base_t, return_t(args_t...)> {
 public:
  /// Base class type
  using base_type = base_t;
  /// Function return type
  using return_type = return_t;
  /// Self type
  using self_type = PairTypeDispatcher<base_type, return_type(args_t...)>;
  /// Internal function signature: both dispatch arguments plus trailing args.
  /// The two leading arguments are always passed in *registration* order.
  using function_signature = return_t(const base_t&, const base_t&, args_t...);
  /// Function pointer type of the trailing (non-dispatched) arguments only
  using function_pointer_type = return_t (*)(args_t...);
  /// Function wrapper type
  using function_type = std::function<function_signature>;

  /// Default constructor
  PairTypeDispatcher() = default;

  /// Constructor that registers multiple function pointers with auto-detected
  /// types.
  template <typename... derived_a_t, typename... derived_b_t>
  explicit PairTypeDispatcher(return_t (*... funcs)(const derived_a_t&,
                                                    const derived_b_t&,
                                                    args_t...))
    requires((std::is_base_of_v<base_t, derived_a_t> && ...) &&
             (std::is_base_of_v<base_t, derived_b_t> && ...))
  {
    (registerFunction(funcs), ...);
  }

  /// Register a free function with explicit or auto-detected derived types.
  /// @tparam derived_a_t The first derived type the function handles
  /// @tparam derived_b_t The second derived type the function handles
  /// @param func The function pointer
  /// @return Reference to this dispatcher for chaining
  template <typename derived_a_t, typename derived_b_t>
    requires std::is_base_of_v<base_t, derived_a_t> &&
             std::is_base_of_v<base_t, derived_b_t>
  self_type& registerFunction(return_t (*func)(const derived_a_t&,
                                               const derived_b_t&, args_t...)) {
    std::type_index typeA(typeid(derived_a_t));
    std::type_index typeB(typeid(derived_b_t));

    // Reject an exact duplicate for this unordered pair
    for (const auto& reg : m_registrations) {
      if (sameUnorderedPair(reg.typeA, reg.typeB, typeA, typeB)) {
        throw std::runtime_error(std::format(
            "Function already registered for type pair: ({}, {})",
            boost::core::demangle(typeA.name()),
            boost::core::demangle(typeB.name())));
      }
    }

    // Best-effort registration-time conflict detection: if both derived types
    // are default constructible we can materialise sample objects and check
    // whether an existing registration would also claim this pair (in either
    // orientation).
    if constexpr (std::is_default_constructible_v<derived_a_t> &&
                  std::is_default_constructible_v<derived_b_t>) {
      derived_a_t sampleA{};
      derived_b_t sampleB{};
      for (const auto& reg : m_registrations) {
        if (reg.matches(sampleA, sampleB).has_value()) {
          throw std::runtime_error(std::format(
              "Registration conflict: type pair ({}, {}) would also be handled "
              "by existing function registered for pair: ({}, {})",
              boost::core::demangle(typeA.name()),
              boost::core::demangle(typeB.name()),
              boost::core::demangle(reg.typeA.name()),
              boost::core::demangle(reg.typeB.name())));
        }
      }
    }

    Registration reg;
    reg.typeA = typeA;
    reg.typeB = typeB;
    reg.checkerA = [](const base_t& obj) {
      return dynamic_cast<const derived_a_t*>(&obj) != nullptr;
    };
    reg.checkerB = [](const base_t& obj) {
      return dynamic_cast<const derived_b_t*>(&obj) != nullptr;
    };
    // The invoker always receives its two leading arguments in registration
    // order: `a` is expected to be a `derived_a_t`, `b` a `derived_b_t`.
    reg.invoke = [func]<typename... Ts>(const base_t& a, const base_t& b,
                                        Ts&&... args) -> return_t {
      const auto* da = dynamic_cast<const derived_a_t*>(&a);
      const auto* db = dynamic_cast<const derived_b_t*>(&b);
      if (da == nullptr || db == nullptr) {
        throw std::bad_cast();
      }
      return func(*da, *db, std::forward<Ts>(args)...);
    };

    m_registrations.push_back(std::move(reg));
    return *this;
  }

  /// Call the registered function for the given pair of objects.
  ///
  /// The pair is matched against every registration in both orientations. The
  /// selected function is always invoked with the arguments normalised to the
  /// order in which it was registered.
  ///
  /// @param lhs The first object to dispatch on
  /// @param rhs The second object to dispatch on
  /// @param args Additional trailing arguments forwarded to the function
  /// @return The return value from the registered function
  template <typename... func_args_t>
  return_t operator()(const base_t& lhs, const base_t& rhs,
                      func_args_t&&... args) const
    requires std::invocable<function_pointer_type, func_args_t...>
  {
    const Registration* match = nullptr;
    bool swap = false;

    for (const auto& reg : m_registrations) {
      auto orientation = reg.matches(lhs, rhs);
      if (!orientation.has_value()) {
        continue;
      }
      if (match != nullptr) {
        throw std::runtime_error(std::format(
            "Ambiguous dispatch for type pair ({}, {}): multiple functions can "
            "handle it",
            boost::core::demangle(typeid(lhs).name()),
            boost::core::demangle(typeid(rhs).name())));
      }
      match = &reg;
      swap = *orientation;
    }

    if (match == nullptr) {
      throw std::runtime_error(
          std::format("No function registered for type pair: ({}, {})",
                      boost::core::demangle(typeid(lhs).name()),
                      boost::core::demangle(typeid(rhs).name())));
    }

    if (swap) {
      return match->invoke(rhs, lhs, std::forward<func_args_t>(args)...);
    }
    return match->invoke(lhs, rhs, std::forward<func_args_t>(args)...);
  }

  /// Check if a function is registered for the given pair of objects.
  /// @param lhs The first object to check
  /// @param rhs The second object to check
  /// @return true if a function is registered, false otherwise
  bool hasFunction(const base_t& lhs, const base_t& rhs) const {
    for (const auto& reg : m_registrations) {
      if (reg.matches(lhs, rhs).has_value()) {
        return true;
      }
    }
    return false;
  }

  /// Check if a function is registered for the given type pair.
  /// @tparam derived_a_t The first type to check
  /// @tparam derived_b_t The second type to check
  /// @return true if a function is registered, false otherwise
  template <typename derived_a_t, typename derived_b_t>
    requires std::is_base_of_v<base_t, derived_a_t> &&
             std::is_base_of_v<base_t, derived_b_t>
  bool hasFunction() const {
    std::type_index typeA(typeid(derived_a_t));
    std::type_index typeB(typeid(derived_b_t));
    for (const auto& reg : m_registrations) {
      if (sameUnorderedPair(reg.typeA, reg.typeB, typeA, typeB)) {
        return true;
      }
    }
    return false;
  }

  /// Clear all registered functions
  void clear() { m_registrations.clear(); }

  /// Get the number of registered functions
  /// @return Number of registered functions
  std::size_t size() const { return m_registrations.size(); }

 private:
  using cast_checker_type = std::function<bool(const base_t&)>;

  struct Registration {
    std::type_index typeA{typeid(void)};
    std::type_index typeB{typeid(void)};
    cast_checker_type checkerA;
    cast_checker_type checkerB;
    function_type invoke;

    /// Determine whether this registration handles the pair `(lhs, rhs)`.
    /// @return std::nullopt if it does not match, otherwise a bool that is
    ///         false for the forward orientation (lhs->A, rhs->B) and true for
    ///         the swapped orientation (rhs->A, lhs->B). The forward
    ///         orientation is preferred when both match.
    std::optional<bool> matches(const base_t& lhs, const base_t& rhs) const {
      if (checkerA(lhs) && checkerB(rhs)) {
        return false;
      }
      if (checkerA(rhs) && checkerB(lhs)) {
        return true;
      }
      return std::nullopt;
    }
  };

  static bool sameUnorderedPair(const std::type_index& a1,
                                const std::type_index& b1,
                                const std::type_index& a2,
                                const std::type_index& b2) {
    return (a1 == a2 && b1 == b2) || (a1 == b2 && b1 == a2);
  }

  std::vector<Registration> m_registrations;
};

}  // namespace Acts

// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Navigation/INavigationPolicy.hpp"
#include "Acts/Navigation/MultiNavigationPolicy.hpp"

#include <concepts>
#include <memory>

namespace Acts {

class TrackingVolume;
class GeometryContext;
class Logger;
class INavigationPolicy;

namespace detail {

/// Trait detecting a `std::unique_ptr<T>` where `T` derives from
/// @ref INavigationPolicy.
/// @tparam T The type to inspect
template <typename T>
struct IsNavigationPolicyPtr : std::false_type {};

template <typename T>
  requires(std::derived_from<T, INavigationPolicy>)
struct IsNavigationPolicyPtr<std::unique_ptr<T>> : std::true_type {};

/// Concept for factory functions that return a concrete navigation policy by
/// value.
/// @tparam F The factory function type
/// @tparam Args The argument types for the factory function
template <typename F, typename... Args>
concept NavigationPolicyValueFactoryConcept = requires(
    F f, const GeometryContext& gctx, const TrackingVolume& volume,
    const Logger& logger, Args&&... args) {
  { f(gctx, volume, logger, args...) } -> std::derived_from<INavigationPolicy>;

  requires NavigationPolicyConcept<decltype(f(gctx, volume, logger, args...))>;
};

/// Concept for factory functions that return a `std::unique_ptr` to a
/// navigation policy. This allows the factory function to select between
/// different concrete policy types at runtime (e.g. based on the volume
/// geometry).
/// @tparam F The factory function type
/// @tparam Args The argument types for the factory function
template <typename F, typename... Args>
concept NavigationPolicyPointerFactoryConcept = requires(
    F f, const GeometryContext& gctx, const TrackingVolume& volume,
    const Logger& logger, Args&&... args) {
  requires IsNavigationPolicyPtr<
      std::remove_cvref_t<decltype(f(gctx, volume, logger, args...))>>::value;
};
}  // namespace detail

/// Base class for navigation policy factories. The factory can be assembled
/// iteratively by using `make` followed by a number of calls to the `add`
/// function of the helper type. Example:
///
/// ```cpp
/// auto factory = NavigationPolicyFactory{}
///  .add<NavigationPolicy1>(arg1, arg2)
///  .add<NavigationPolicy2>(/*no args*/)
///  .asUniquePtr();
/// ```
class NavigationPolicyFactory {
 private:
  /// Type alias for factory functions that create navigation policies
  using factory_type = std::function<std::unique_ptr<INavigationPolicy>(
      const GeometryContext&, const TrackingVolume&, const Logger&)>;

  /// Private constructor for internal use
  /// @param factories Vector of factory functions
  explicit NavigationPolicyFactory(std::vector<factory_type>&& factories)
      : m_factories(std::move(factories)) {}

 public:
  /// Default constructor
  NavigationPolicyFactory() = default;

  /// Add a navigation policy to the factory
  /// @tparam P The policy type to add
  /// @param args The arguments to pass to the policy constructor
  /// @note Arguments need to be copy constructible because the factory must be
  ///       able to execute multiple times.
  /// @return New instance of this object with the added factory for method
  ///         chaining
  template <NavigationPolicyConcept P, typename... Args>
    requires(std::is_constructible_v<P, const GeometryContext&,
                                     const TrackingVolume&, const Logger&,
                                     Args...> &&
             (std::is_copy_constructible_v<Args> && ...))
  constexpr NavigationPolicyFactory add(Args&&... args) && {
    auto factory = [=](const GeometryContext& gctx,
                       const TrackingVolume& volume, const Logger& logger) {
      return std::make_unique<P>(gctx, volume, logger, args...);
    };

    m_factories.push_back(std::move(factory));
    return std::move(*this);
  }

  /// Add a policy created by a factory function that returns a
  /// `std::unique_ptr` to a navigation policy. Returning a pointer allows the
  /// factory function to select between different concrete policy types at
  /// runtime (e.g. based on the volume geometry).
  /// @tparam Fn The type of the function to construct the policy
  /// @param fn The factory function
  /// @param args The arguments to pass to the policy factory
  /// @note Arguments need to be copy constructible because the factory must be
  ///       able to execute multiple times.
  /// @return New instance of this object with the added factory for method
  ///         chaining
  template <typename Fn, typename... Args>
    requires(detail::NavigationPolicyPointerFactoryConcept<Fn, Args...> &&
             (std::is_copy_constructible_v<Args> && ...))
  constexpr NavigationPolicyFactory add(Fn&& fn, Args&&... args) && {
    auto factory = [=](const GeometryContext& gctx,
                       const TrackingVolume& volume,
                       const Logger& logger) -> factory_type::result_type {
      return fn(gctx, volume, logger, args...);
    };

    m_factories.push_back(std::move(factory));
    return std::move(*this);
  }

  /// Add a policy created by a factory function that returns a concrete policy
  /// by value.
  /// @tparam Fn The type of the function to construct the policy
  /// @param fn The factory function
  /// @param args The arguments to pass to the policy factory
  /// @note Arguments need to be copy constructible because the factory must be
  ///       able to execute multiple times.
  /// @return New instance of this object with the added factory for method
  ///         chaining
  /// @deprecated Return a `std::unique_ptr<INavigationPolicy>` from the factory
  ///             function instead. This also lets the function select the
  ///             concrete policy type at runtime.
  template <typename Fn, typename... Args>
    requires(detail::NavigationPolicyValueFactoryConcept<Fn, Args...> &&
             !detail::NavigationPolicyPointerFactoryConcept<Fn, Args...> &&
             (std::is_copy_constructible_v<Args> && ...))
  [[deprecated(
      "Return a std::unique_ptr<INavigationPolicy> from the factory function "
      "instead of a policy by value.")]]
  constexpr NavigationPolicyFactory add(Fn&& fn, Args&&... args) && {
    auto factory = [=](const GeometryContext& gctx,
                       const TrackingVolume& volume,
                       const Logger& logger) -> factory_type::result_type {
      using policy_type = decltype(fn(gctx, volume, logger, args...));
      return std::make_unique<policy_type>(fn(gctx, volume, logger, args...));
    };

    m_factories.push_back(std::move(factory));
    return std::move(*this);
  }

  /// Move the factory into a unique pointer
  /// @return A unique pointer to the factory
  std::unique_ptr<NavigationPolicyFactory> asUniquePtr() && {
    return std::make_unique<NavigationPolicyFactory>(std::move(*this));
  }

  /// Construct a multi-navigation policy using the registered factories
  /// @param gctx The geometry context
  /// @param volume The tracking volume
  /// @param logger The logger
  /// @return A unique pointer to the constructed MultiNavigationPolicy
  /// @throws std::runtime_error if no factories are registered
  std::unique_ptr<MultiNavigationPolicy> operator()(
      const GeometryContext& gctx, const TrackingVolume& volume,
      const Logger& logger) const {
    if (m_factories.empty()) {
      throw std::runtime_error(
          "No factories registered in the navigation policy factory");
    }

    std::vector<std::unique_ptr<INavigationPolicy>> policies;
    policies.reserve(m_factories.size());
    for (auto& factory : m_factories) {
      policies.push_back(factory(gctx, volume, logger));
    }

    return std::make_unique<MultiNavigationPolicy>(std::move(policies));
  }

  /// Construct a navigation policy using the factories (alias for operator())
  /// @param gctx The geometry context
  /// @param volume The tracking volume
  /// @param logger The logger
  /// @return A unique pointer to the constructed navigation policy
  std::unique_ptr<INavigationPolicy> build(const GeometryContext& gctx,
                                           const TrackingVolume& volume,
                                           const Logger& logger) const {
    return operator()(gctx, volume, logger);
  }

 private:
  /// Vector of factory functions to create navigation policies
  std::vector<factory_type> m_factories;
};

}  // namespace Acts

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
#include <functional>
#include <memory>
#include <span>
#include <stdexcept>
#include <vector>

namespace Acts {

class TrackingVolume;
class GeometryContext;
class Logger;
class INavigationPolicy;

/// Non-owning view of the policies already constructed by the same
/// @ref NavigationPolicyFactory, in the order they were added. Passed to
/// policies that build on an earlier sibling. The pointees are owned by the
/// @ref MultiNavigationPolicy the factory returns, so they outlive the policy
/// receiving them; the view itself is only valid during construction.
using NavigationPolicySiblings = std::span<const INavigationPolicy* const>;

namespace detail {

/// Trait detecting a `std::unique_ptr<T>` where `T` derives from
/// @ref INavigationPolicy.
/// @tparam T The type to inspect
template <typename T>
struct IsNavigationPolicyPtr : std::false_type {};

template <typename T>
  requires(std::derived_from<T, INavigationPolicy>)
struct IsNavigationPolicyPtr<std::unique_ptr<T>> : std::true_type {};

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

/// Same as @ref NavigationPolicyPointerFactoryConcept, for factory functions
/// that additionally receive the earlier sibling policies.
/// @tparam F The factory function type
/// @tparam Args The argument types for the factory function
template <typename F, typename... Args>
concept NavigationPolicySiblingPointerFactoryConcept = requires(
    F f, const GeometryContext& gctx, const TrackingVolume& volume,
    const Logger& logger, NavigationPolicySiblings siblings, Args&&... args) {
  requires IsNavigationPolicyPtr<std::remove_cvref_t<decltype(f(
      gctx, volume, logger, siblings, args...))>>::value;
};

/// Policy constructible from the volume and @p Args
template <typename P, typename... Args>
concept NavigationPolicyConstructible =
    std::is_constructible_v<P, const GeometryContext&, const TrackingVolume&,
                            const Logger&, Args...>;

/// Policy constructible from the volume, the earlier siblings and @p Args
template <typename P, typename... Args>
concept NavigationPolicySiblingConstructible =
    std::is_constructible_v<P, const GeometryContext&, const TrackingVolume&,
                            const Logger&, NavigationPolicySiblings, Args...>;
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
///
/// Policies are constructed in the order they were added. A policy (or factory
/// function) may take the already constructed earlier policies as a
/// @ref NavigationPolicySiblings argument right after the logger, which lets it
/// build on an earlier sibling. Siblings are owned by the same
/// @ref MultiNavigationPolicy, so references to them stay valid.
class NavigationPolicyFactory {
 private:
  /// Type alias for factory functions that create navigation policies
  using factory_type = std::function<std::unique_ptr<INavigationPolicy>(
      const GeometryContext&, const TrackingVolume&, const Logger&,
      NavigationPolicySiblings)>;

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
  /// @note The policy is constructed from `(gctx, volume, logger, args...)`,
  ///       or, if it accepts them, from
  ///       `(gctx, volume, logger, siblings, args...)` with the policies
  ///       constructed before it.
  /// @note Arguments need to be copy constructible because the factory must be
  ///       able to execute multiple times.
  /// @return New instance of this object with the added factory for method
  ///         chaining
  template <NavigationPolicyConcept P, typename... Args>
    requires((detail::NavigationPolicyConstructible<P, Args...> ||
              detail::NavigationPolicySiblingConstructible<P, Args...>) &&
             (std::is_copy_constructible_v<Args> && ...))
  constexpr NavigationPolicyFactory add(Args&&... args) && {
    auto factory =
        [=](const GeometryContext& gctx, const TrackingVolume& volume,
            const Logger& logger,
            NavigationPolicySiblings siblings) -> factory_type::result_type {
      // Prefer the sibling-aware constructor: a policy offering it wants them
      if constexpr (detail::NavigationPolicySiblingConstructible<P, Args...>) {
        return std::make_unique<P>(gctx, volume, logger, siblings, args...);
      } else {
        return std::make_unique<P>(gctx, volume, logger, args...);
      }
    };

    m_factories.push_back(std::move(factory));
    return std::move(*this);
  }

  /// Add a policy created by a factory function that returns a
  /// `std::unique_ptr` to a navigation policy. Returning a pointer allows the
  /// factory function to select between different concrete policy types at
  /// runtime (e.g. based on the volume geometry).
  /// @tparam Fn The type of the function to construct the policy
  /// @param fn The factory function, called with
  ///           `(gctx, volume, logger, args...)` or, if it accepts them,
  ///           `(gctx, volume, logger, siblings, args...)`
  /// @param args The arguments to pass to the policy factory
  /// @note Arguments need to be copy constructible because the factory must be
  ///       able to execute multiple times.
  /// @return New instance of this object with the added factory for method
  ///         chaining
  template <typename Fn, typename... Args>
    requires(
        (detail::NavigationPolicyPointerFactoryConcept<Fn, Args...> ||
         detail::NavigationPolicySiblingPointerFactoryConcept<Fn, Args...>) &&
        (std::is_copy_constructible_v<Args> && ...))
  constexpr NavigationPolicyFactory add(Fn&& fn, Args&&... args) && {
    auto factory =
        [=](const GeometryContext& gctx, const TrackingVolume& volume,
            const Logger& logger,
            NavigationPolicySiblings siblings) -> factory_type::result_type {
      if constexpr (detail::NavigationPolicySiblingPointerFactoryConcept<
                        Fn, Args...>) {
        return fn(gctx, volume, logger, siblings, args...);
      } else {
        return fn(gctx, volume, logger, args...);
      }
    };

    m_factories.push_back(std::move(factory));
    return std::move(*this);
  }

  /// Move the factory into a unique pointer
  /// @return A unique pointer to the factory
  std::unique_ptr<NavigationPolicyFactory> asUniquePtr() && {
    return std::make_unique<NavigationPolicyFactory>(std::move(*this));
  }

  /// Construct a multi-navigation policy using the registered factories. The
  /// factories run in order, each receiving the policies built before it.
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
    // Observation-only view handed to the children, hides the ownership
    std::vector<const INavigationPolicy*> siblings;
    siblings.reserve(m_factories.size());
    for (const auto& factory : m_factories) {
      policies.push_back(
          factory(gctx, volume, logger, NavigationPolicySiblings{siblings}));
      siblings.push_back(policies.back().get());
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

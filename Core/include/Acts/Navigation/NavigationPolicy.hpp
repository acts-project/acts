// This file is part of the Acts project.
//
// Copyright (C) 2024 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Navigation/NavigationDelegate.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/DelegateChain.hpp"

#include <tuple>
#include <type_traits>
#include <utility>

namespace Acts {

class TrackingVolume;
class INavigationPolicy;

template <typename T>
concept NavigationPolicyConcept = requires {
  requires std::is_base_of_v<INavigationPolicy, T>;
  // Has a conforming update method
  requires requires(T policy, const NavigationArguments& args) {
    { policy.updateState(args) };
  };
};

class INavigationPolicy {
 public:
  static void noopUpdate(const NavigationArguments& /*unused*/) {}

  virtual ~INavigationPolicy() = default;

  virtual void connect(NavigationDelegate& delegate) const = 0;

 protected:
  template <NavigationPolicyConcept T>
  void connectDefault(NavigationDelegate& delegate) const {
    // This cannot be a concept because we use it in CRTP below
    const auto* self = static_cast<const T*>(this);
    DelegateChainFactory{delegate}.add<&T::updateState>(self).store(delegate);
  }
};

class TryAllPortalNavigationPolicy final : public INavigationPolicy {
 public:
  explicit TryAllPortalNavigationPolicy(const TrackingVolume& volume)
      : m_volume(&volume) {}

  // @NOTE: This implementation is PRELIMINARY! It is subject to change!
  void updateState(const NavigationArguments& args) const;

  void connect(NavigationDelegate& delegate) const override {
    connectDefault<TryAllPortalNavigationPolicy>(delegate);
  }

 private:
  const TrackingVolume* m_volume;
};

static_assert(NavigationPolicyConcept<TryAllPortalNavigationPolicy>);

class MultiNavigationPolicyBase : public INavigationPolicy {};

template <NavigationPolicyConcept... Policies>
class MultiNavigationPolicy final : public MultiNavigationPolicyBase {
 public:
  MultiNavigationPolicy(Policies&&... policies)
      : m_policies{std::forward<Policies>(policies)...} {}

  void connect(NavigationDelegate& delegate) const override {
    auto factory = add(DelegateChainFactory{delegate},
                       std::index_sequence_for<Policies...>{});

    factory.store(delegate);
  }

  const std::tuple<Policies...>& policies() const { return m_policies; }

 private:
  template <std::size_t... Is>
  constexpr auto add(
      auto factory,
      std::integer_sequence<std::size_t, Is...> /*unused*/) const {
    return add<Is...>(factory);
  }

  template <std::size_t I, std::size_t... Is>
  constexpr auto add(auto factory) const {
    using policy_type = std::tuple_element_t<I, decltype(m_policies)>;
    auto* policy = &std::get<I>(m_policies);

    auto added = factory.template add<&policy_type::updateState>(policy);

    if constexpr (sizeof...(Is) > 0) {
      return add<Is...>(added);
    } else {
      return added;
    }
  }

  std::tuple<Policies...> m_policies;
};

template <typename... Factories>
class NavigationPolicyFactory {
  template <typename... Fs>
  using policy_type =
      std::conditional_t<(sizeof...(Fs) > 0),
                         decltype(MultiNavigationPolicy(std::invoke(
                             std::declval<Fs>(),
                             std::declval<const TrackingVolume&>())...)),
                         MultiNavigationPolicy<>>;

 public:
  template <typename...>
  friend class NavigationPolicyFactory;

  NavigationPolicyFactory() = default;

  template <NavigationPolicyConcept P, typename... Args>
  constexpr auto add(Args&&... args)
    requires(std::is_constructible_v<P, const TrackingVolume&, Args...>)
  {
    auto factory = [&](const TrackingVolume& volume) {
      return P{volume, std::forward<Args>(args)...};
    };

    return NavigationPolicyFactory<Factories..., decltype(factory)>{
        std::tuple_cat(std::move(m_factories),
                       std::make_tuple(std::move(factory)))};
  }

  auto operator()(const TrackingVolume& volume) const
    requires(sizeof...(Factories) > 0)
  {
    using policy_type = decltype(MultiNavigationPolicy(std::invoke(
        std::declval<Factories>(), std::declval<const TrackingVolume&>())...));

    return std::apply(
        [&](auto&&... factories) -> std::unique_ptr<policy_type> {
          return std::make_unique<policy_type>(
              std::invoke(factories, volume)...);
        },
        m_factories);
  }

 private:
  NavigationPolicyFactory(std::tuple<Factories...>&& factories)
      : m_factories(std::move(factories)) {}

  std::tuple<Factories...> m_factories;
};

}  // namespace Acts

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

#include <memory>
namespace Acts {

class TrackingVolume;
class INavigationPolicy;

namespace detail {
template <typename... Factories>
class NavigationPolicyFactoryImpl;
}

class NavigationPolicyFactory {
 public:
  virtual ~NavigationPolicyFactory() = default;

  // This needs to be listed here, but the return type cannot be spelled out
  // yet.
  static auto make();

  // This will potentially get serialization interface and deserialization
  // functionality

  virtual std::unique_ptr<INavigationPolicy> build(
      const TrackingVolume& volume) const = 0;
};

namespace detail {

template <>
class NavigationPolicyFactoryImpl<> {
 public:
  template <typename...>
  friend class NavigationPolicyFactoryImpl;
  NavigationPolicyFactoryImpl() = default;

  // Arguments need to be copy constructible because the factory must be able to
  // execute multiple times.
  template <NavigationPolicyConcept P, typename... Args>
  constexpr auto add(Args&&... args)
    requires(std::is_constructible_v<P, const TrackingVolume&, Args...> &&
             (std::is_copy_constructible_v<Args> && ...))
  {
    auto factory = [&](const TrackingVolume& volume) {
      return P{volume, std::forward<Args>(args)...};
    };

    return NavigationPolicyFactoryImpl<decltype(factory)>{
        std::tuple_cat(std::make_tuple(std::move(factory)))};
  }
};

template <typename F, typename... Fs>
class NavigationPolicyFactoryImpl<F, Fs...> : public NavigationPolicyFactory {
  template <typename... Ts>
  using policy_type_helper = decltype(MultiNavigationPolicy(std::invoke(
      std::declval<Ts>(), std::declval<const TrackingVolume&>())...));

  using policy_type = policy_type_helper<F, Fs...>;

 public:
  template <NavigationPolicyConcept P, typename... Args>
  constexpr auto add(Args&&... args)
    requires(std::is_constructible_v<P, const TrackingVolume&, Args...> &&
             (std::is_copy_constructible_v<Args> && ...))
  {
    auto factory = [&](const TrackingVolume& volume) {
      return P{volume, std::forward<Args>(args)...};
    };

    return NavigationPolicyFactoryImpl<F, Fs..., decltype(factory)>{
        std::tuple_cat(std::move(m_factories),
                       std::make_tuple(std::move(factory)))};
  }

  constexpr std::unique_ptr<NavigationPolicyFactoryImpl<F, Fs...>>
  asUniquePtr() && {
    return std::make_unique<NavigationPolicyFactoryImpl<F, Fs...>>(
        std::move(*this));
  }

  std::unique_ptr<policy_type> operator()(const TrackingVolume& volume) const {
    return std::apply(
        [&](auto&&... factories) -> std::unique_ptr<policy_type> {
          return std::make_unique<policy_type>(
              std::invoke(factories, volume)...);
        },
        m_factories);
  }

  std::unique_ptr<INavigationPolicy> build(
      const TrackingVolume& volume) const override {
    return operator()(volume);
  }

 private:
  template <typename...>
  friend class NavigationPolicyFactoryImpl;

  NavigationPolicyFactoryImpl(std::tuple<F, Fs...>&& factories)
      : m_factories(std::move(factories)) {}

  std::tuple<F, Fs...> m_factories;
};

}  // namespace detail

inline auto NavigationPolicyFactory::make() {
  return detail::NavigationPolicyFactoryImpl<>{};
}

}  // namespace Acts

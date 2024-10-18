// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Navigation/INavigationPolicy.hpp"

namespace Acts {

class MultiNavigationPolicyBase : public INavigationPolicy {};

template <NavigationPolicyConcept... Policies>
  requires(sizeof...(Policies) > 0)
class MultiNavigationPolicy final : public MultiNavigationPolicyBase {
 public:
  MultiNavigationPolicy(Policies&&... policies)
      : m_policies{std::forward<Policies>(policies)...} {}

  void connect(NavigationDelegate& delegate) const override {
    auto factory = add(DelegateChainBuilder{delegate},
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
}  // namespace Acts

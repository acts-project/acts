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

/// Combined navigation policy that calls all contained other navigation
/// policies.
class MultiNavigationPolicy final : public INavigationPolicy {
 public:
  template <typename... Policies>
  explicit MultiNavigationPolicy(std::unique_ptr<Policies>... policies)
      : MultiNavigationPolicy{[](auto... args) {
          std::vector<std::unique_ptr<INavigationPolicy>> policyPtrs;
          auto fill = [&policyPtrs](auto& policy) {
            policyPtrs.push_back(std::move(policy));
          };

          (fill(args), ...);

          return policyPtrs;
        }(std::move(policies)...)} {}

  explicit MultiNavigationPolicy(
      std::vector<std::unique_ptr<INavigationPolicy>>&& policies)
      : m_policyPtrs(std::move(policies)) {
    m_delegates.reserve(m_policyPtrs.size());
    for (auto& policy : m_policyPtrs) {
      auto& delegate = m_delegates.emplace_back();
      policy->connect(delegate);
      if (!delegate.connected()) {
        throw std::runtime_error(
            "Failed to connect policy in MultiNavigationPolicyDynamic");
      }
    }
  }

  void connect(NavigationDelegate& delegate) const override {
    if (!std::ranges::all_of(m_delegates,
                             [](const auto& d) { return d.connected(); })) {
      throw std::runtime_error(
          "Not all delegates are connected to MultiNavigationPolicyDynamic");
    }
    delegate.connect<&MultiNavigationPolicy::initializeCandidates>(this);
  }

  const std::span<const std::unique_ptr<INavigationPolicy>> policies() const {
    return std::span{m_policyPtrs};
  }

 private:
  void initializeCandidates(const NavigationArguments& args,
                            AppendOnlyNavigationStream& stream,
                            const Logger& logger) const {
    for (const auto& delegate : m_delegates) {
      delegate(args, stream, logger);
    }
  }

  std::vector<std::unique_ptr<INavigationPolicy>> m_policyPtrs;

  std::vector<NavigationDelegate> m_delegates;
};
}  // namespace Acts

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
      std::vector<std::unique_ptr<INavigationPolicy>>&& policies);

  void connect(NavigationDelegate& delegate) const override;

  const std::span<const std::unique_ptr<INavigationPolicy>> policies() const;

 private:
  void initializeCandidates(const NavigationArguments& args,
                            AppendOnlyNavigationStream& stream,
                            const Logger& logger) const;

  std::vector<std::unique_ptr<INavigationPolicy>> m_policyPtrs;

  std::vector<NavigationDelegate> m_delegates;
};
}  // namespace Acts

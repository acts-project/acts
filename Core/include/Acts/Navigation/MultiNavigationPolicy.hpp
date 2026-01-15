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
/// policies. This class manages multiple navigation policies and executes
/// them sequentially to populate the navigation stream with candidates.
class MultiNavigationPolicy final : public INavigationPolicy {
 public:
  /// Constructor from multiple unique_ptr policies
  /// @tparam Policies The types of the navigation policies
  /// @param policies Unique pointers to navigation policies
  template <typename... Policies>
  explicit MultiNavigationPolicy(std::unique_ptr<Policies>... policies)
      : MultiNavigationPolicy{[](auto... args) {
          std::vector<std::unique_ptr<INavigationPolicy>> policyPtrs;

          if constexpr (sizeof...(args) > 0) {
            auto fill = [&policyPtrs](auto& policy) {
              policyPtrs.push_back(std::move(policy));
            };

            (fill(args), ...);
          }

          return policyPtrs;
        }(std::move(policies)...)} {}

  /// Constructor from a vector of navigation policies
  /// @param policies Vector of unique pointers to navigation policies
  explicit MultiNavigationPolicy(
      std::vector<std::unique_ptr<INavigationPolicy>>&& policies);

  /// Connect this policy to a navigation delegate
  /// @param delegate The navigation delegate to connect to
  void connect(NavigationDelegate& delegate) const override;

  /// Access the contained navigation policies
  /// @return Span of const unique pointers to the navigation policies
  std::span<const std::unique_ptr<INavigationPolicy>> policies() const;

  /// Visit this policy first, and then all child policies one by one.
  /// @param visitor The function to call for each policy
  void visit(const std::function<void(const INavigationPolicy&)>& visitor)
      const override;

 private:
  /// Initialize navigation candidates by calling all contained policies
  /// @param gctx The geometry context
  /// @param args The navigation arguments
  /// @param stream The navigation stream to populate
  /// @param logger Logger for debug output
  void initializeCandidates(const GeometryContext& gctx,
                            const NavigationArguments& args,
                            NavigationPolicyState& state,
                            AppendOnlyNavigationStream& stream,
                            const Logger& logger) const;

  /// Vector of unique pointers to the contained navigation policies
  std::vector<std::unique_ptr<INavigationPolicy>> m_policyPtrs;

  /// Vector of navigation delegates, one for each policy
  std::vector<NavigationDelegate> m_delegates;
};
}  // namespace Acts

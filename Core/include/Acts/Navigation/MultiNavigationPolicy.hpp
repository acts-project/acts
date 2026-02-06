// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Navigation/INavigationPolicy.hpp"
#include "Acts/Utilities/Zip.hpp"

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

  /// State structure for MultiNavigationPolicy
  /// Holds the states for all contained child policies
  struct State {
    /// Vector of navigation policy states, one for each child policy
    std::vector<NavigationPolicyState> policyStates;
  };

  /// Check if all child policies are in a valid state
  /// @param gctx The geometry context
  /// @param args The navigation arguments
  /// @param state The navigation policy state to check
  /// @param logger Logger for debug output
  /// @return True if all child policy states are valid, false otherwise
  bool isValid(const GeometryContext& gctx, const NavigationArguments args,
               NavigationPolicyState& state,
               const Logger& logger) const override {
    ACTS_VERBOSE("MultiNavigationPolicy isValid check, forward to "
                 << m_policyPtrs.size() << " policies.");

    auto& thisState = state.as<State>();
    if (thisState.policyStates.size() != m_policyPtrs.size()) {
      ACTS_ERROR("MultiNavigationPolicy isValid: number of states ("
                 << thisState.policyStates.size()
                 << ") does not match number of policies ("
                 << m_policyPtrs.size() << ").");
      throw std::runtime_error(
          "MultiNavigationPolicy isValid: inconsistent state size.");
    }

    for (auto [policy, policyState] :
         zip(m_policyPtrs, thisState.policyStates)) {
      if (!policy->isValid(gctx, args, policyState, logger)) {
        return false;
      }
    }
    // If no policy rejected, return true
    return true;
  }

  /// Create and initialize states for this policy and all child policies
  /// @param gctx The geometry context
  /// @param args The navigation arguments
  /// @param stateManager The state manager to push the new states onto
  /// @param logger Logger for debug output
  void createState(const GeometryContext& gctx, const NavigationArguments args,
                   NavigationPolicyStateManager& stateManager,
                   const Logger& logger) const override {
    ACTS_VERBOSE("MultiNavigationPolicy createState, create states for "
                 << m_policyPtrs.size() << " policies.");

    // Push child states first, then at the end push this policy's state,
    // containing the references

    std::vector<NavigationPolicyState> states;
    states.reserve(m_policyPtrs.size());

    for (const auto& policy : m_policyPtrs) {
      ACTS_VERBOSE("Creating child state for policy ");
      policy->createState(gctx, args, stateManager, logger);
      states.emplace_back(stateManager.currentState());
    }

    ACTS_VERBOSE("Created "
                 << states.size()
                 << " child states for MultiNavigationPolicy (of which "
                 << std::ranges::count_if(
                        states, [](const auto& s) { return !s.empty(); })
                 << " are non-empty)");

    // Important, push at the end
    stateManager.pushState<State>(std::move(states));
  }

  /// Remove the states for this policy and all child policies from the state manager
  /// @param stateManager The state manager to pop the states from
  /// @param logger Logger for debug output
  void popState(NavigationPolicyStateManager& stateManager,
                const Logger& logger) const override {
    // By default, we didn't push anything, so we don't need to pop anything
    ACTS_VERBOSE("MultiNavigationPolicy popState called, popping for "
                 << m_policyPtrs.size() << " child policies");

    // Pops this policy's state first
    stateManager.popState();

    // Then pops all child states
    for (const auto& policy : m_policyPtrs) {
      policy->popState(stateManager, logger);
    }
  }

 private:
  /// Initialize navigation candidates by calling all contained policies
  /// @param gctx The geometry context
  /// @param args The navigation arguments
  /// @param state The navigation policy state
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

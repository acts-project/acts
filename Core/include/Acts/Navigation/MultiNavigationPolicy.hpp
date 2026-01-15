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

#include <ranges>

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

  struct State {
    std::vector<NavigationPolicyState> policyStates;
  };

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

  NavigationPolicyState createState(const GeometryContext& gctx,
                                    const NavigationArguments args,
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
      states.emplace_back(
          policy->createState(gctx, args, stateManager, logger));
    }

    ACTS_VERBOSE("Created "
                 << states.size()
                 << " child states for MultiNavigationPolicy (of which "
                 << std::ranges::count_if(
                        states, [](const auto& s) { return !s.empty(); })
                 << " are valid)");

    auto [state, any] = stateManager.pushState<State>(std::move(states));
    return any;
  }

  void popState(NavigationPolicyStateManager& stateManager,
                const Logger& logger) const override {
    // By default, we didn't push anything, so we don't need to poop anything
    ACTS_VERBOSE("MultiNavigationPolicy popState called, popping for "
                 << m_policyPtrs.size() << " child policies");

    // `createState` pushed a State containing all child states, so we
    // pop the current state (which should be this policy's state) and then pop
    // however many non-empty states the children had pushed.

    std::vector<NavigationPolicyState> states =
        std::move(stateManager.currentState().as<State>().policyStates);
    stateManager.popState();

    if (states.size() != m_policyPtrs.size()) {
      ACTS_ERROR("MultiNavigationPolicy popState: number of states ("
                 << states.size() << ") does not match number of policies ("
                 << m_policyPtrs.size() << ").");
      throw std::runtime_error(
          "MultiNavigationPolicy popState: inconsistent state size.");
    }

    // @TODO: Possibly use updated zip | reverse
    for (std::size_t i = states.size(); i-- > 0;) {
      auto& state = states[i];
      auto& policy = m_policyPtrs[i];
      if (!state.empty()) {
        policy->popState(stateManager, logger);
      }
    }
  }

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

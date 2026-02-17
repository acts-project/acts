// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Navigation/MultiNavigationPolicy.hpp"

namespace Acts {

MultiNavigationPolicy::MultiNavigationPolicy(
    std::vector<std::unique_ptr<INavigationPolicy>>&& policies)
    : m_policyPtrs(std::move(policies)) {
  m_delegates.reserve(m_policyPtrs.size());
  for (const auto& policy : m_policyPtrs) {
    auto& delegate = m_delegates.emplace_back();
    policy->connect(delegate);
    if (!delegate.connected()) {
      throw std::runtime_error(
          "Failed to connect policy in MultiNavigationPolicyDynamic");
    }
  }
}

void MultiNavigationPolicy::connect(NavigationDelegate& delegate) const {
  if (!std::ranges::all_of(m_delegates,
                           [](const auto& d) { return d.connected(); })) {
    throw std::runtime_error(
        "Not all delegates are connected to MultiNavigationPolicyDynamic");
  }
  delegate.connect<&MultiNavigationPolicy::initializeCandidates>(this);
}

std::span<const std::unique_ptr<INavigationPolicy>>
MultiNavigationPolicy::policies() const {
  return std::span{m_policyPtrs};
}

void MultiNavigationPolicy::initializeCandidates(
    [[maybe_unused]] const GeometryContext& gctx,
    const NavigationArguments& args, NavigationPolicyState& state,
    AppendOnlyNavigationStream& stream, const Logger& logger) const {
  auto& thisState = state.as<State>();
  if (thisState.policyStates.size() != m_delegates.size()) {
    ACTS_ERROR("MultiNavigationPolicy initializeCandidates: number of states ("
               << thisState.policyStates.size()
               << ") does not match number of policies (" << m_delegates.size()
               << ").");
    throw std::runtime_error(
        "MultiNavigationPolicy initializeCandidates: inconsistent state size.");
  }

  for (auto [delegate, policyState] :
       zip(m_delegates, thisState.policyStates)) {
    delegate(gctx, args, policyState, stream, logger);
  }
}

void MultiNavigationPolicy::visit(
    const std::function<void(const INavigationPolicy&)>& visitor) const {
  visitor(*this);
  for (const auto& policy : m_policyPtrs) {
    visitor(*policy);
  }
}

void MultiNavigationPolicy::createState(
    const GeometryContext& gctx, const NavigationArguments& args,
    NavigationPolicyStateManager& stateManager, const Logger& logger) const {
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

  ACTS_VERBOSE(
      "Created " << states.size()
                 << " child states for MultiNavigationPolicy (of which "
                 << std::ranges::count_if(
                        states, [](const auto& s) { return !s.empty(); })
                 << " are non-empty)");

  // Important, push at the end
  stateManager.pushState<State>(std::move(states));
}

void MultiNavigationPolicy::popState(NavigationPolicyStateManager& stateManager,
                                     const Logger& logger) const {
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

bool MultiNavigationPolicy::isValid(const GeometryContext& gctx,
                                    const NavigationArguments& args,
                                    NavigationPolicyState& state,
                                    const Logger& logger) const {
  ACTS_VERBOSE("MultiNavigationPolicy isValid check, forward to "
               << m_policyPtrs.size() << " policies.");

  auto& thisState = state.as<State>();
  if (thisState.policyStates.size() != m_policyPtrs.size()) {
    ACTS_ERROR("MultiNavigationPolicy isValid: number of states ("
               << thisState.policyStates.size()
               << ") does not match number of policies (" << m_policyPtrs.size()
               << ").");
    throw std::runtime_error(
        "MultiNavigationPolicy isValid: inconsistent state size.");
  }

  for (auto [policy, policyState] : zip(m_policyPtrs, thisState.policyStates)) {
    if (!policy->isValid(gctx, args, policyState, logger)) {
      return false;
    }
  }
  // If no policy rejected, return true
  return true;
}

}  // namespace Acts

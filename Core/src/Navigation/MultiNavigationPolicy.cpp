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
  // Child states were pushed contiguously below this policy's own state (see
  // createState); their count always equals the number of policies. When the
  // navigator skipped state creation for a stateless volume, the state handle
  // is the empty sentinel; the (equally stateless) children then receive
  // empty handles too, which they by definition do not read.
  const bool haveStates = !state.empty();
  const std::size_t childBase =
      haveStates ? state.index() - m_delegates.size() : 0;
  for (std::size_t i = 0; i < m_delegates.size(); ++i) {
    NavigationPolicyState childState =
        haveStates ? state.atIndex(childBase + i) : NavigationPolicyState{};
    m_delegates[i](gctx, args, childState, stream, logger);
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

  // Push each child state onto the manager. They end up contiguously on the
  // stack and are recovered by index. Track whether any child pushed a real
  // (non-default) state, which is the signal that this composite can invalidate.
  bool anyChildConstrained = false;
  for (const auto& policy : m_policyPtrs) {
    ACTS_VERBOSE("Creating child state for policy ");
    policy->createState(gctx, args, stateManager, logger);
    if (!stateManager.currentState().isDefault()) {
      anyChildConstrained = true;
    }
  }

  // Push this policy's own state on top of its children. Use the default
  // EmptyState -- the "no validity constraint" sentinel -- unless some child
  // can invalidate, in which case push the marker State so the navigator
  // consults isValid() for this volume.
  if (anyChildConstrained) {
    stateManager.pushState<State>();
  } else {
    stateManager.pushState<EmptyState>();
  }
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

  // Child states were pushed contiguously below this policy's own state (see
  // createState); their count always equals the number of policies.
  const std::size_t childBase = state.index() - m_policyPtrs.size();
  for (std::size_t i = 0; i < m_policyPtrs.size(); ++i) {
    NavigationPolicyState childState = state.atIndex(childBase + i);
    if (!m_policyPtrs[i]->isValid(gctx, args, childState, logger)) {
      return false;
    }
  }
  // If no policy rejected, return true
  return true;
}

}  // namespace Acts

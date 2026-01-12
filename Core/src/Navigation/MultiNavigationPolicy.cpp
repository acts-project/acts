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
    const NavigationArguments& args, AppendOnlyNavigationStream& stream,
    const Logger& logger) const {
  for (const auto& delegate : m_delegates) {
    delegate(gctx, args, stream, logger);
  }
}

void MultiNavigationPolicy::visit(
    const std::function<void(const INavigationPolicy&)>& visitor) const {
  visitor(*this);
  for (const auto& policy : m_policyPtrs) {
    visitor(*policy);
  }
}

}  // namespace Acts

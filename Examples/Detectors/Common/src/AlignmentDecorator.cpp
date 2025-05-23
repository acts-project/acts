// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/DetectorCommons/AlignmentDecorator.hpp"

#include "Acts/Geometry/AlignmentDelegate.hpp"
#include "Acts/Geometry/TransformStore.hpp"

#include <ranges>

ActsExamples::AlignmentDecorator::AlignmentDecorator(
    const Config& cfg, std::unique_ptr<const Acts::Logger> logger)
    : m_cfg(cfg), m_logger(std::move(logger)) {
  if (m_cfg.alignmentStores.empty() && m_cfg.nominalStore == nullptr) {
    throw std::invalid_argument(
        "Missing alignment stores (and nominal store), run without alignment "
        "decorator!");
  }
  // Sort on leading IOV
  std::ranges::sort(m_cfg.alignmentStores,
                    [](const auto& lhs, const auto& rhs) {
                      const auto& [lhsIov, lhsStore] = lhs;
                      const auto& [rhsIov, rhsStore] = rhs;
                      return lhsIov[0u] < rhsIov[0u];
                    });
  // Check for overlapping IOVs
  for (const auto [istore, iovStore] : Acts::enumerate(m_cfg.alignmentStores)) {
    if (istore > 0) {
      const auto& [iov, store] = iovStore;
      const auto& [prevIov, prevStore] = m_cfg.alignmentStores[istore - 1];
      if (iov[0] == prevIov[0] || prevIov[1] >= iov[0]) {
        throw std::invalid_argument(
            "Intersecting IOVs found as [" + std::to_string(prevIov[0]) + ", " +
            std::to_string(prevIov[1]) + "] and [" + std::to_string(iov[0]) +
            ", " + std::to_string(iov[1]) + "]");
      }
    }
  }
}

ActsExamples::ProcessCode ActsExamples::AlignmentDecorator::decorate(
    AlgorithmContext& context) {
  std::size_t eventNumber = context.eventNumber;

  // Start with the current alignment store
  auto currentStore = m_cfg.nominalStore;
  auto matchedStore = std::ranges::find_if(
      m_cfg.alignmentStores, [eventNumber](const auto& iovStore) {
        const auto& [iov, store] = iovStore;
        return iov[0] >= eventNumber && eventNumber <= iov[1];
      });
  if (matchedStore != m_cfg.alignmentStores.end()) {
    const auto& [iov, store] = *matchedStore;
    ACTS_VERBOSE("Found alignment store for event number " +
                 std::to_string(eventNumber) + " in [" +
                 std::to_string(iov[0]) + ", " + std::to_string(iov[1]) + "]");
    currentStore = store;
  }

  // We must have a valid alignment store at this point
  if (currentStore == nullptr) {
    throw std::invalid_argument(
        "No alignment store found for event number " +
        std::to_string(eventNumber) +
        ", check IOV bounds and/or configuration of nominal alignment store");
  }

  // Create a DetectorElement alignment store for this context
  Acts::AlignmentDelegate currentAlignment;
  currentAlignment.connect<&Acts::ITransformStore::contextualTransform>(
      currentStore.get());
  context.geoContext = Acts::GeometryContext{currentAlignment};

  return ActsExamples::ProcessCode::SUCCESS;
}

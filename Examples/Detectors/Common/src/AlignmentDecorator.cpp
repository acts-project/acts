// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/DetectorCommons/AlignmentDecorator.hpp"

#include "ActsExamples/DetectorCommons/AlignmentContext.hpp"

#include <ranges>
#include <stdexcept>

namespace ActsExamples {

namespace {

template <typename IovCollection>
void checkOverlappingIOVs(IovCollection& iovStores) {
  // First sort
  std::ranges::sort(iovStores, [](const auto& a, const auto& b) {
    return std::get<0>(a)[0] < std::get<0>(b)[0];
  });
  // Then check for overlaps
  for (std::size_t i = 1; i < iovStores.size(); ++i) {
    const auto& prevIov = std::get<0>(iovStores[i - 1]);
    const auto& currIov = std::get<0>(iovStores[i]);
    if (currIov[0] < prevIov[1]) {
      throw std::runtime_error(
          "Overlapping IOVs detected: " + std::to_string(prevIov[0]) + "-" +
          std::to_string(prevIov[1]) + " and " + std::to_string(currIov[0]) +
          "-" + std::to_string(currIov[1]));
    }
  }
}

bool eventWithinIOV(const std::array<std::size_t, 2>& iov,
                    std::size_t eventNumber) {
  return eventNumber >= iov[0] && eventNumber <= iov[1];
}

}  // namespace

AlignmentDecorator::AlignmentDecorator(const AlignmentDecorator::Config& config,
                                       Acts::Logging::Level level)
    : m_cfg(config),
      m_logger(Acts::getDefaultLogger("AlignmentDecorator", level)) {
  // Check for overlapping IOVs in the stores
  if (!m_cfg.iovStores.empty() && m_cfg.iovGenerators.empty()) {
    ACTS_VERBOSE("Configuring AlignmentDecorator with IOV stores.");
    checkOverlappingIOVs(m_cfg.iovStores);
    m_iovStores.reserve(m_cfg.iovStores.size());
    std::ranges::for_each(m_cfg.iovStores, [&](const auto& iovStore) {
      const auto& [iov, store] = iovStore;
      m_iovStores.emplace_back(iov, store, 0);  // Initialize counter to 0
    });
  } else if (!m_cfg.iovGenerators.empty()) {
    ACTS_VERBOSE("Configuring AlignmentDecorator with IOV generators.");
    m_iovGenerators = m_cfg.iovGenerators;
    checkOverlappingIOVs(m_cfg.iovGenerators);
  } else {
    throw std::runtime_error(
        "Both IOV stores and generators cannot be either empty nor "
        "configured.");
  }
}

ProcessCode AlignmentDecorator::decorate(AlgorithmContext& context) {
  // Get the event number & look for the alignment store
  const std::size_t eventNumber = context.eventNumber;
  // A pointer to an alignment store - if set it will be used to decorate the
  // geometry context consequently
  IAlignmentStore* alignmentStore = nullptr;
  if (!m_iovStores.empty()) {
    ACTS_VERBOSE("Looking for alignment store for event " << eventNumber
                                                          << " in IOV stores.");
    std::scoped_lock lock(m_iovMutex);
    std::ranges::for_each(m_iovStores, [&](auto& iovStore) {
      auto& [iov, store, counter] = iovStore;
      if (eventWithinIOV(iov, eventNumber)) {
        alignmentStore = store.get();
        if (counter == 0u) {
          ++counter;
        }  // Increase the counter if first use;
      }
    });
  }

  if (!m_iovGenerators.empty() && alignmentStore == nullptr) {
    ACTS_VERBOSE("No alignment store found for event "
                 << eventNumber << ", checking generators.");
    std::ranges::for_each(m_iovGenerators, [&](auto& iovGenerator) {
      const auto& [iov, generator] = iovGenerator;
      if (eventWithinIOV(iov, eventNumber)) {
        auto alignmentStorePtr = m_cfg.nominalStore->clone();
        std::scoped_lock lock(m_iovMutex);
        alignmentStorePtr->visitStore(generator);
        m_iovStores.emplace_back(iov, alignmentStorePtr,
                                 1);  // Add to the IOV stores
        alignmentStore = alignmentStorePtr.get();
        ACTS_DEBUG("Generated alignment store for event "
                   << eventNumber << " with IOV [" << iov[0] << ", " << iov[1]
                   << "]");
      }
    });
  }

  // Run the garbage collection if requested
  if (m_cfg.garbageCollection) {
    std::scoped_lock lock(m_iovMutex);
    // First increase the counters
    for (auto& [iov, store, counter] : m_iovStores) {
      if (!eventWithinIOV(iov, eventNumber) && counter > 0) {
        ++counter;
        ACTS_VERBOSE("Garbage collection: keeping alignment store for IOV ["
                     << iov[0] << ", " << iov[1]
                     << "] events passed since last in use: " << counter);
      }
    }
    // Remove if count is greater than the interval
    auto rmIdx = std::ranges::remove_if(m_iovStores, [this](
                                                         const auto& iovStore) {
      const auto& [iov, store, counter] = iovStore;
      if (counter > m_cfg.gcInterval) {
        ACTS_DEBUG("Garbage collection: removing alignment store for IOV [" +
                   std::to_string(iov[0]) + ", " + std::to_string(iov[1]) +
                   "] with counter " + std::to_string(counter - 1) +
                   " reaching/exceeding interval " +
                   std::to_string(m_cfg.gcInterval));
        return true;  // Remove if counter exceeds interval
      }
      return false;  // Keep it
    });
    // Do the actual removal
    m_iovStores.erase(rmIdx.begin(), rmIdx.end());
  }

  // Decorate the alignment context if a store is set
  if (alignmentStore != nullptr) {
    ACTS_DEBUG("Decorating AlgorithmContext with alignment store for event "
               << eventNumber);
    AlignmentContext alignmentContext{alignmentStore};
    context.geoContext = Acts::GeometryContext(alignmentContext);
  } else {
    ACTS_VERBOSE("No alignment store found for event "
                 << eventNumber << ", skipping decoration.");
  }

  return ProcessCode::SUCCESS;
}

}  // namespace ActsExamples

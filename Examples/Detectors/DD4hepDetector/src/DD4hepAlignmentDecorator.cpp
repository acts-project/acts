// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/DD4hepDetector/DD4hepAlignmentDecorator.hpp"

#include "Acts/Plugins/DD4hep/DD4hepAlignmentStore.hpp"
#include "Acts/Plugins/DD4hep/DD4hepDetectorElement.hpp"

ActsExamples::DD4hepAlignmentDecorator::DD4hepAlignmentDecorator(
    const Config& cfg, std::unique_ptr<const Acts::Logger> logger)
    : m_cfg(cfg), m_logger(std::move(logger)) {
  if (m_cfg.alignmentStores.empty() && m_cfg.nominalStore == nullptr) {
    throw std::invalid_argument(
        "Missing alignment stores (and nominal store), run without alignment "
        "decorator!");
  }
}

ActsExamples::ProcessCode ActsExamples::DD4hepAlignmentDecorator::decorate(
    AlgorithmContext& context) {
  std::size_t eventNumber = context.eventNumber;

  // Start with the current alignment store
  auto currentStore = m_cfg.nominalStore;
  if (eventNumber > m_cfg.alignmentStores.begin()->first) {
    auto currentStoreBound = m_cfg.alignmentStores.lower_bound(eventNumber);
    if (currentStoreBound != m_cfg.alignmentStores.end()) {
      currentStore = currentStoreBound->second;
    }
  }
  // We must have a valid alignment store at this point
  if (currentStore == nullptr) {
    throw std::invalid_argument(
        "No alignment store found for event number " +
        std::to_string(eventNumber) +
        ", check IOV bounds and/or configuration of nominal alignment store");
  }

  // Create a DetectorElement alignment store for this context
  Acts::DD4hepGeometryContext::Alignment currentAlignment; 
  currentAlignment.connect<&Acts::IDD4hepAlignmentStore::contextualTransform>(currentStore.get());
  context.geoContext = Acts::DD4hepGeometryContext(currentAlignment);
  return ActsExamples::ProcessCode::SUCCESS;
}

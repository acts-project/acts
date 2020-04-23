// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ACTFW/ContextualDetector/PayloadDecorator.hpp"

#include "ACTFW/ContextualDetector/PayloadDetectorElement.hpp"
#include "Acts/Geometry/TrackingGeometry.hpp"

FW::Contextual::PayloadDecorator::PayloadDecorator(
    const FW::Contextual::PayloadDecorator::Config& cfg,
    std::unique_ptr<const Acts::Logger> logger)
    : m_cfg(cfg), m_logger(std::move(logger)) {
  if (m_cfg.trackingGeometry != nullptr) {
    // parse and populate
    parseGeometry(*m_cfg.trackingGeometry.get());
  }
}

FW::ProcessCode FW::Contextual::PayloadDecorator::decorate(
    AlgorithmContext& context) {
  // Start with the nominal map
  std::vector<Acts::Transform3D> aStore = m_nominalStore;

  ACTS_VERBOSE("New IOV detected, emulate new alignment");
  if (context.eventNumber % m_cfg.iovSize) {
    for (auto& tf : aStore) {
      tf *= Acts::AngleAxis3D(m_cfg.rotationStep * context.eventNumber,
                              Acts::Vector3D::UnitY());
    }
    // This creates a full payload context, i.e. the nominal store
    PayloadDetectorElement::ContextType alignableGeoContext;
    alignableGeoContext.alignmentStore = std::move(aStore);
    context.geoContext =
        std::make_any<PayloadDetectorElement::ContextType>(alignableGeoContext);
  }
  return ProcessCode::SUCCESS;
}

void FW::Contextual::PayloadDecorator::parseGeometry(
    const Acts::TrackingGeometry& tGeometry) {
  // Double-visit - first count
  size_t nTransforms = 0;
  tGeometry.visitSurfaces([&nTransforms](const auto*) { ++nTransforms; });

  PayloadDetectorElement::ContextType nominalCtx;

  // Collect the surfacas into the nominal store
  std::vector<Acts::Transform3D> aStore(nTransforms,
                                        Acts::Transform3D::Identity());

  auto fillTransforms = [&aStore, &nominalCtx](const auto* surface) -> void {
    auto alignableElement = dynamic_cast<const PayloadDetectorElement*>(
        surface->associatedDetectorElement());
    aStore[alignableElement->identifier()] = surface->transform(nominalCtx);
  };

  tGeometry.visitSurfaces(fillTransforms);
  m_nominalStore = std::move(aStore);
}

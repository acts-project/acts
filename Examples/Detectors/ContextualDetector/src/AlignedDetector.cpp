// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/ContextualDetector/AlignedDetector.hpp"

#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/ContextualDetector/AlignmentDecorator.hpp"
#include "ActsExamples/ContextualDetector/ExternalAlignmentDecorator.hpp"
#include "ActsExamples/ContextualDetector/ExternallyAlignedDetectorElement.hpp"
#include "ActsExamples/ContextualDetector/InternalAlignmentDecorator.hpp"
#include "ActsExamples/ContextualDetector/InternallyAlignedDetectorElement.hpp"
#include "ActsExamples/Framework/RandomNumbers.hpp"

namespace ActsExamples {

AlignedDetector::AlignedDetector(const Config& cfg)
    : GenericDetector(cfg, NoBuildTag{}), m_cfg(cfg) {
  if (m_cfg.mode == Config::Mode::External) {
    ExternallyAlignedDetectorElement::ContextType nominalContext;
    m_nominalGeometryContext = Acts::GeometryContext(nominalContext);
  } else {
    InternallyAlignedDetectorElement::ContextType nominalContext;
    nominalContext.nominal = true;
    m_nominalGeometryContext = Acts::GeometryContext(nominalContext);
  }

  // Let's create a random number service
  RandomNumbers::Config randomNumberConfig;
  randomNumberConfig.seed = m_cfg.seed;
  auto randomNumberSvc = std::make_shared<RandomNumbers>(randomNumberConfig);

  auto fillDecoratorConfig = [&](AlignmentDecorator::Config& config) {
    config.iovSize = m_cfg.iovSize;
    config.flushSize = m_cfg.flushSize;
    config.doGarbageCollection = m_cfg.doGarbageCollection;

    // The misalignments
    config.gSigmaX = m_cfg.sigmaInPlane;
    config.gSigmaY = m_cfg.sigmaInPlane;
    config.gSigmaZ = m_cfg.sigmaOutPlane;
    config.aSigmaX = m_cfg.sigmaOutRot;
    config.aSigmaY = m_cfg.sigmaOutRot;
    config.aSigmaZ = m_cfg.sigmaInRot;
    config.randomNumberSvc = randomNumberSvc;
    config.firstIovNominal = m_cfg.firstIovNominal;
  };

  if (m_cfg.mode == Config::Mode::External) {
    ExternalAlignmentDecorator::Config agcsConfig;
    fillDecoratorConfig(agcsConfig);

    auto detectorElementFactory =
        [this](std::shared_ptr<const Acts::Transform3> transform,
               std::shared_ptr<const Acts::PlanarBounds> bounds,
               double thickness,
               std::shared_ptr<const Acts::ISurfaceMaterial> material)
        -> std::shared_ptr<GenericDetectorElement> {
      auto id = m_detectorStore.size();
      auto detElem = std::make_shared<ExternallyAlignedDetectorElement>(
          id, std::move(transform), std::move(bounds), thickness,
          std::move(material));
      m_detectorStore.push_back(detElem);
      return detElem;
    };

    buildTrackingGeometry(detectorElementFactory);

    agcsConfig.trackingGeometry = m_trackingGeometry;

    m_contextDecorators.push_back(std::make_shared<ExternalAlignmentDecorator>(
        std::move(agcsConfig),
        Acts::getDefaultLogger("AlignmentDecorator", m_cfg.decoratorLogLevel)));
  } else {
    InternalAlignmentDecorator::Config agcsConfig;
    fillDecoratorConfig(agcsConfig);

    auto detectorElementFactory =
        [this, &agcsConfig](
            std::shared_ptr<const Acts::Transform3> transform,
            std::shared_ptr<const Acts::PlanarBounds> bounds, double thickness,
            std::shared_ptr<const Acts::ISurfaceMaterial> material)
        -> std::shared_ptr<GenericDetectorElement> {
      auto id = m_detectorStore.size();
      auto detElem = std::make_shared<InternallyAlignedDetectorElement>(
          id, std::move(transform), std::move(bounds), thickness,
          std::move(material));
      m_detectorStore.push_back(detElem);
      agcsConfig.detectorStore.push_back(detElem);
      return detElem;
    };

    buildTrackingGeometry(detectorElementFactory);

    m_contextDecorators.push_back(std::make_shared<InternalAlignmentDecorator>(
        std::move(agcsConfig),
        Acts::getDefaultLogger("AlignmentDecorator", m_cfg.decoratorLogLevel)));
  }
}

}  // namespace ActsExamples

// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/ContextualDetector/AlignedDetector.hpp"

#include "Acts/Geometry/DetectorElementBase.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/ContextualDetector/AlignmentDecorator.hpp"
#include "ActsExamples/ContextualDetector/ExternalAlignmentDecorator.hpp"
#include "ActsExamples/ContextualDetector/ExternallyAlignedDetectorElement.hpp"
#include "ActsExamples/ContextualDetector/InternalAlignmentDecorator.hpp"
#include "ActsExamples/ContextualDetector/InternallyAlignedDetectorElement.hpp"
#include "ActsExamples/DetectorCommons/DetectorBase.hpp"
#include "ActsExamples/Framework/RandomNumbers.hpp"
#include "ActsExamples/GenericDetector/BuildGenericDetector.hpp"

#include <memory>

namespace ActsExamples {

AlignedDetectorFactory::AlignedDetectorFactory(const Config& cfg)
    : DetectorFactoryBase(
          Acts::getDefaultLogger("AlignedDetectorFactory", cfg.logLevel)),
      m_cfg(cfg) {}

std::shared_ptr<DetectorBase> AlignedDetectorFactory::buildDetector() const {
  Acts::GeometryContext geometryContext;
  std::vector<std::shared_ptr<const Acts::DetectorElementBase>> detectorStore;
  std::shared_ptr<const Acts::TrackingGeometry> gen1Geometry;
  std::shared_ptr<Acts::Experimental::Detector> gen2Geometry;
  std::vector<std::shared_ptr<IContextDecorator>> contextDecorators;

  if (m_cfg.mode == Config::Mode::External) {
    InternallyAlignedDetectorElement::ContextType nominalContext;
    geometryContext = Acts::GeometryContext(nominalContext);
  } else {
    InternallyAlignedDetectorElement::ContextType nominalContext;
    nominalContext.nominal = true;
    geometryContext = Acts::GeometryContext(nominalContext);
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

    std::vector<std::vector<std::shared_ptr<ExternallyAlignedDetectorElement>>>
        specificDetectorStore;

    gen1Geometry = Generic::buildDetector<ExternallyAlignedDetectorElement>(
        geometryContext, specificDetectorStore, m_cfg.buildLevel,
        m_cfg.materialDecorator, m_cfg.buildProto, m_cfg.surfaceLogLevel,
        m_cfg.layerLogLevel, m_cfg.volumeLogLevel);
    agcsConfig.trackingGeometry = gen1Geometry;

    // need to upcast to store in this object as well
    for (auto& lstore : specificDetectorStore) {
      for (auto& ldet : lstore) {
        detectorStore.push_back(ldet);
      }
    }

    contextDecorators.push_back(std::make_shared<ExternalAlignmentDecorator>(
        std::move(agcsConfig),
        Acts::getDefaultLogger("AlignmentDecorator", m_cfg.decoratorLogLevel)));
  } else {
    InternalAlignmentDecorator::Config agcsConfig;
    fillDecoratorConfig(agcsConfig);

    gen1Geometry = Generic::buildDetector<InternallyAlignedDetectorElement>(
        geometryContext, agcsConfig.detectorStore, m_cfg.buildLevel,
        m_cfg.materialDecorator, m_cfg.buildProto, m_cfg.surfaceLogLevel,
        m_cfg.layerLogLevel, m_cfg.volumeLogLevel);

    // need to upcast to store in this object as well
    for (auto& lstore : agcsConfig.detectorStore) {
      for (auto& ldet : lstore) {
        detectorStore.push_back(ldet);
      }
    }

    contextDecorators.push_back(std::make_shared<InternalAlignmentDecorator>(
        std::move(agcsConfig),
        Acts::getDefaultLogger("AlignmentDecorator", m_cfg.decoratorLogLevel)));
  }

  return std::make_shared<PreConstructedDetector>(
      std::move(geometryContext), std::move(detectorStore),
      std::move(gen1Geometry), std::move(gen2Geometry),
      std::move(contextDecorators));
}

}  // namespace ActsExamples

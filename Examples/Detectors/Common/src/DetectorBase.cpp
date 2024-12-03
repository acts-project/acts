// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/DetectorCommons/DetectorBase.hpp"

#include "Acts/Detector/Detector.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/TrackingGeometry.hpp"
#include "ActsExamples/Framework/IContextDecorator.hpp"

namespace ActsExamples {

DetectorBase::DetectorBase() = default;

DetectorBase::DetectorBase(
    std::vector<std::shared_ptr<const Acts::DetectorElementBase>> detectorStore,
    std::vector<std::shared_ptr<IContextDecorator>> contextDecorators)
    : m_detectorStore(std::move(detectorStore)),
      m_contextDecorators(std::move(contextDecorators)) {}

DetectorBase::~DetectorBase() = default;

std::vector<std::shared_ptr<IContextDecorator>>
DetectorBase::contextDecorators() const {
  return m_contextDecorators;
}

std::unique_ptr<G4VUserDetectorConstruction>
DetectorBase::buildGeant4DetectorConstruction(
    const Geant4ConstructionOptions& /*options*/) const {
  return nullptr;
}

PreConstructedDetector::PreConstructedDetector() = default;

PreConstructedDetector::PreConstructedDetector(
    Acts::GeometryContext geometryContext,
    std::vector<std::shared_ptr<const Acts::DetectorElementBase>> detectorStore,
    std::shared_ptr<const Acts::TrackingGeometry> trackingGeometry,
    std::shared_ptr<Acts::Experimental::Detector> gen2Geometry,
    std::vector<std::shared_ptr<IContextDecorator>> contextDecorators)
    : DetectorBase(std::move(detectorStore), std::move(contextDecorators)),
      m_geometryContext(std::move(geometryContext)),
      m_trackingGeometry(std::move(trackingGeometry)),
      m_gen2Geometry(std::move(gen2Geometry)) {}

PreConstructedDetector::~PreConstructedDetector() = default;

const Acts::GeometryContext& PreConstructedDetector::geometryContext() const {
  return m_geometryContext;
}

std::shared_ptr<const Acts::TrackingGeometry>
PreConstructedDetector::trackingGeometry() const {
  return m_trackingGeometry;
}

std::shared_ptr<Acts::Experimental::Detector>
PreConstructedDetector::gen2Geometry() const {
  return m_gen2Geometry;
}

DetectorFactoryBase::DetectorFactoryBase(
    std::unique_ptr<const Acts::Logger> logger)
    : m_logger(std::move(logger)) {}

DetectorFactoryBase::DetectorFactoryBase(DetectorFactoryBase&&) = default;

DetectorFactoryBase::~DetectorFactoryBase() = default;

DetectorFactoryBase& DetectorFactoryBase::operator=(DetectorFactoryBase&&) =
    default;

const Acts::Logger& DetectorFactoryBase::logger() const {
  return *m_logger;
}

}  // namespace ActsExamples

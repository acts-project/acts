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
#include "ActsExamples/DetectorCommons/Geant4DetectorConstructionFactory.hpp"
#include "ActsExamples/Framework/IContextDecorator.hpp"

namespace ActsExamples {

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

PreConstructedDetector::PreConstructedDetector(
    Acts::GeometryContext geometryContext,
    std::vector<std::shared_ptr<const Acts::DetectorElementBase>> detectorStore,
    std::shared_ptr<const Acts::TrackingGeometry> gen1Geometry,
    std::shared_ptr<Acts::Experimental::Detector> gen2Geometry,
    std::vector<std::shared_ptr<IContextDecorator>> contextDecorators,
    std::shared_ptr<Geant4DetectorConstructionFactory>
        geant4DetectorConstructionFactory)
    : m_geometryContext(std::move(geometryContext)),
      m_detectorStore(std::move(detectorStore)),
      m_gen1Geometry(std::move(gen1Geometry)),
      m_gen2Geometry(std::move(gen2Geometry)),
      m_contextDecorators(std::move(contextDecorators)),
      m_geant4DetectorConstructionFactory(
          std::move(geant4DetectorConstructionFactory)) {}

PreConstructedDetector::~PreConstructedDetector() = default;

const Acts::GeometryContext& PreConstructedDetector::geometryContext() const {
  return m_geometryContext;
}

const std::vector<std::shared_ptr<const Acts::DetectorElementBase>>&
PreConstructedDetector::detectorStore() const {
  return m_detectorStore;
}

const std::shared_ptr<const Acts::TrackingGeometry>&
PreConstructedDetector::gen1Geometry() const {
  return m_gen1Geometry;
}

const std::shared_ptr<Acts::Experimental::Detector>&
PreConstructedDetector::gen2Geometry() const {
  return m_gen2Geometry;
}

const std::vector<std::shared_ptr<IContextDecorator>>&
PreConstructedDetector::contextDecorators() const {
  return m_contextDecorators;
}

const std::shared_ptr<Geant4DetectorConstructionFactory>&
PreConstructedDetector::geant4DetectorConstructionFactory() const {
  return m_geant4DetectorConstructionFactory;
}

}  // namespace ActsExamples

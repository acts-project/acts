// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Geometry/GeometryContext.hpp"

#include <memory>
#include <vector>

namespace Acts {
class GeometryContext;
class TrackingGeometry;
class DetectorElementBase;
class IMaterialDecorator;
class Logger;
namespace Experimental {
class Detector;
}  // namespace Experimental
}  // namespace Acts

namespace ActsExamples {
class IContextDecorator;
class Geant4DetectorConstructionFactory;
namespace Geant4 {
class RegionCreator;
}  // namespace Geant4
}  // namespace ActsExamples

namespace ActsExamples {

class DetectorBase {
 public:
  virtual ~DetectorBase() = default;

  virtual const Acts::GeometryContext& geometryContext() const = 0;

  virtual const std::vector<std::shared_ptr<const Acts::DetectorElementBase>>&
  detectorStore() const = 0;
  virtual const std::shared_ptr<const Acts::TrackingGeometry>& gen1Geometry()
      const = 0;
  virtual const std::shared_ptr<Acts::Experimental::Detector>& gen2Geometry()
      const = 0;
  virtual const std::vector<std::shared_ptr<ActsExamples::IContextDecorator>>&
  contextDecorators() const = 0;

  virtual const std::shared_ptr<Geant4DetectorConstructionFactory>&
  geant4DetectorConstructionFactory() const = 0;
};

class DetectorFactoryBase {
 public:
  explicit DetectorFactoryBase(std::unique_ptr<const Acts::Logger> logger);
  DetectorFactoryBase(const DetectorFactoryBase&) = delete;
  DetectorFactoryBase(DetectorFactoryBase&&);
  virtual ~DetectorFactoryBase();
  DetectorFactoryBase& operator=(const DetectorFactoryBase&) = delete;
  DetectorFactoryBase& operator=(DetectorFactoryBase&&);

  virtual std::shared_ptr<DetectorBase> buildDetector() const = 0;

 protected:
  std::unique_ptr<const Acts::Logger> m_logger;

  const Acts::Logger& logger() const;
};

class PreConstructedDetector : public DetectorBase {
 public:
  PreConstructedDetector(
      Acts::GeometryContext geometryContext,
      std::vector<std::shared_ptr<const Acts::DetectorElementBase>>
          detectorStore,
      std::shared_ptr<const Acts::TrackingGeometry> gen1Geometry,
      std::shared_ptr<Acts::Experimental::Detector> gen2Geometry,
      std::vector<std::shared_ptr<ActsExamples::IContextDecorator>>
          contextDecorators,
      std::shared_ptr<Geant4DetectorConstructionFactory>
          geant4DetectorConstructionFactory);
  ~PreConstructedDetector() override;

  const Acts::GeometryContext& geometryContext() const override;

  const std::vector<std::shared_ptr<const Acts::DetectorElementBase>>&
  detectorStore() const override;

  const std::shared_ptr<const Acts::TrackingGeometry>& gen1Geometry()
      const override;

  const std::shared_ptr<Acts::Experimental::Detector>& gen2Geometry()
      const override;

  const std::vector<std::shared_ptr<ActsExamples::IContextDecorator>>&
  contextDecorators() const override;

  const std::shared_ptr<Geant4DetectorConstructionFactory>&
  geant4DetectorConstructionFactory() const override;

 private:
  Acts::GeometryContext m_geometryContext;
  std::vector<std::shared_ptr<const Acts::DetectorElementBase>> m_detectorStore;
  std::shared_ptr<const Acts::TrackingGeometry> m_gen1Geometry;
  std::shared_ptr<Acts::Experimental::Detector> m_gen2Geometry;
  std::vector<std::shared_ptr<ActsExamples::IContextDecorator>>
      m_contextDecorators;
  std::shared_ptr<Geant4DetectorConstructionFactory>
      m_geant4DetectorConstructionFactory;
};

}  // namespace ActsExamples

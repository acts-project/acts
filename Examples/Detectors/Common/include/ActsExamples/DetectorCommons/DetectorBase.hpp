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
struct Geant4ConstructionOptions;
}  // namespace ActsExamples

class G4VUserDetectorConstruction;

namespace ActsExamples {

/// Base class for detector instances
class DetectorBase {
 public:
  DetectorBase();
  DetectorBase(
      std::vector<std::shared_ptr<const Acts::DetectorElementBase>>
          detectorStore,
      std::vector<std::shared_ptr<IContextDecorator>> contextDecorators);
  virtual ~DetectorBase();

  virtual const Acts::GeometryContext& geometryContext() const = 0;

  virtual std::shared_ptr<const Acts::TrackingGeometry> trackingGeometry()
      const = 0;
  virtual std::shared_ptr<Acts::Experimental::Detector> gen2Geometry()
      const = 0;
  virtual std::vector<std::shared_ptr<IContextDecorator>> contextDecorators()
      const;

  /// Build the Geant4 detector construction
  /// @note This throws an exception if Geant4 is not enabled
  /// @param options The Geant4 construction options
  /// @return The Geant4 detector construction
  virtual std::unique_ptr<G4VUserDetectorConstruction>
  buildGeant4DetectorConstruction(
      const Geant4ConstructionOptions& options) const;

 protected:
  std::vector<std::shared_ptr<const Acts::DetectorElementBase>> m_detectorStore;
  std::vector<std::shared_ptr<IContextDecorator>> m_contextDecorators;
};

/// Base class for detector factories
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
  PreConstructedDetector();
  PreConstructedDetector(
      Acts::GeometryContext geometryContext,
      std::vector<std::shared_ptr<const Acts::DetectorElementBase>>
          detectorStore,
      std::shared_ptr<const Acts::TrackingGeometry> trackingGeometry,
      std::shared_ptr<Acts::Experimental::Detector> gen2Geometry,
      std::vector<std::shared_ptr<IContextDecorator>> contextDecorators);
  ~PreConstructedDetector() override;

  const Acts::GeometryContext& geometryContext() const override;

  std::shared_ptr<const Acts::TrackingGeometry> trackingGeometry()
      const override;

  std::shared_ptr<Acts::Experimental::Detector> gen2Geometry() const override;

 private:
  Acts::GeometryContext m_geometryContext;
  std::shared_ptr<const Acts::TrackingGeometry> m_trackingGeometry;
  std::shared_ptr<Acts::Experimental::Detector> m_gen2Geometry;
};

}  // namespace ActsExamples

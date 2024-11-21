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

class G4VUserDetectorConstruction;

namespace ActsExamples {

struct Gen1GeometryHolder {
  Acts::GeometryContext geometryContext;

  std::shared_ptr<const Acts::TrackingGeometry> trackingGeometry;
  std::vector<std::shared_ptr<ActsExamples::IContextDecorator>>
      contextDecorators;
  std::vector<std::shared_ptr<const Acts::DetectorElementBase>> detectorStore;
};

struct Gen2GeometryHolder {
  Acts::GeometryContext geometryContext;

  std::shared_ptr<Acts::Experimental::Detector> detector;
  std::vector<std::shared_ptr<ActsExamples::IContextDecorator>>
      contextDecorators;
  std::vector<std::shared_ptr<const Acts::DetectorElementBase>> detectorStore;
};

class DetectorBase {
 public:
  explicit DetectorBase(std::unique_ptr<const Acts::Logger> logger);
  DetectorBase(const DetectorBase&) = delete;
  DetectorBase(DetectorBase&&);
  virtual ~DetectorBase();
  DetectorBase& operator=(const DetectorBase&) = delete;
  DetectorBase& operator=(DetectorBase&&);

  virtual Gen1GeometryHolder buildGen1Geometry();

  virtual Gen2GeometryHolder buildGen2Geometry();

  virtual std::shared_ptr<Geant4DetectorConstructionFactory>
  buildGeant4DetectorConstructionFactory();

 protected:
  std::unique_ptr<const Acts::Logger> m_logger;

  const Acts::Logger& logger() const { return *m_logger; }
};

}  // namespace ActsExamples

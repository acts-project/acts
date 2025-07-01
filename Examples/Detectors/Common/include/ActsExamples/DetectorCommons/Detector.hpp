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

///  @brief Baseline class to represent a complete experimental setup. The Detector
///         is constructed by an external geometry plugin (e.g. DD4Hep,
///         GeoModel) and holds the tracking geometry representation & the full
///         Geant4 geometry description.
class Detector {
 public:
  /// @brief Standard constructor taking a logger object
  explicit Detector(std::unique_ptr<const Acts::Logger> logger);
  virtual ~Detector();
  /// @brief Returns the reference to the empty unaligned geometry context
  virtual const Acts::GeometryContext& nominalGeometryContext() const;
  /// @brief Returns the valid pointer to the tracking geometry.
  ///        Throws an exception if the pointer is invalid
  virtual std::shared_ptr<const Acts::TrackingGeometry> trackingGeometry()
      const;
  /// @brief Returns the valid pointer to the Gen-2 tracking geometry.
  ///        Throws an exception if the pointer is invalid
  virtual std::shared_ptr<Acts::Experimental::Detector> gen2Geometry() const;
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
  const Acts::Logger& logger() const;

  std::unique_ptr<const Acts::Logger> m_logger;

  Acts::GeometryContext m_nominalGeometryContext;
  std::shared_ptr<const Acts::TrackingGeometry> m_trackingGeometry;
  std::shared_ptr<Acts::Experimental::Detector> m_gen2Geometry;
  std::vector<std::shared_ptr<const Acts::DetectorElementBase>> m_detectorStore;
  std::vector<std::shared_ptr<IContextDecorator>> m_contextDecorators;
};

}  // namespace ActsExamples

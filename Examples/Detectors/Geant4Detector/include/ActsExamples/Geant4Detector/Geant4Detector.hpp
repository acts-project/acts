// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Detector/ProtoDetector.hpp"
#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "Acts/Plugins/Geant4/Geant4DetectorSurfaceFactory.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/DetectorCommons/Detector.hpp"

#include <memory>
#include <string>
#include <tuple>
#include <vector>

class G4VPhysicalVolume;

namespace Acts {
class TrackingGeometry;
class Geant4DetectorElement;
class Surface;

namespace Experimental {
class Detector;
}
}  // namespace Acts

namespace ActsExamples {
class IContextDecorator;
}  // namespace ActsExamples

namespace ActsExamples::Geant4 {

struct Geant4Detector : public DetectorCommons::Detector {
  using ContextDecorators =
      std::vector<std::shared_ptr<ActsExamples::IContextDecorator>>;

  using DetectorElements =
      std::vector<std::shared_ptr<Acts::Geant4DetectorElement>>;
  using Surfaces = std::vector<std::shared_ptr<Acts::Surface>>;

  /// Nested configuration struct
  struct Config {
    /// The detector/geometry name
    std::string name = "";
    /// The Geant4 world volume
    const G4VPhysicalVolume* g4World = nullptr;
    /// The Converter options: detector surfaces
    Acts::Geant4DetectorSurfaceFactory::Options g4SurfaceOptions;
    /// The corresponding ProtoDetector
    Acts::ProtoDetector protoDetector;
    /// Optional geometry identifier hook to be used during closure
    std::shared_ptr<const Acts::GeometryIdentifierHook> geometryIdentifierHook =
        std::make_shared<Acts::GeometryIdentifierHook>();
    /// Logging level of the child tools
    Acts::Logging::Level logLevel = Acts::Logging::INFO;
  };

  explicit Geant4Detector(const Config& cfg);

  const DetectorElements& detectorElements() const;

  void drop() final;

 private:
  Config m_cfg;

  DetectorElements m_detectorElements;

  std::unique_ptr<const Acts::Logger> m_logger;

  const Acts::Logger& logger() const { return *m_logger; }

  void buildTrackingGeometry(const Acts::GeometryContext& gctx) final;
  void buildDetector() final;

  /// @brief Construct an Acts::Detector from a Geant4 world volume
  /// @param cfg the configuration of the Geant4 detector
  /// @param logger a logger instance
  /// @return a tuple of an Acts::Detector object, a ContextDecorator & the created detector elements
  std::tuple<std::shared_ptr<Acts::Experimental::Detector>, ContextDecorators,
             DetectorElements>
  constructDetector() const;

  /// @brief Construct a TrackingGeometry from a Geant4 world volume using the KDTreeTrackingGeometryBuilder builder
  ///
  /// @param cfg the configuration of the Geant4 detector
  /// @param kdtCfg the configuration of the KDTreeTrackingGeometryBuilder
  /// @param logger a logger instance
  ///
  /// @return a tuple of an Acts::TrackingGeometry object,  a ContextDecorator & the created detector elements
  std::tuple<std::shared_ptr<const Acts::TrackingGeometry>, ContextDecorators,
             DetectorElements>
  constructTrackingGeometry() const;

  /// @brief Convert Geant4VPhysicalVolume objects into Acts components
  ///
  /// @param cfg the configuration of the Geant4 detector
  /// @param logger a logger instance
  ///
  /// @return a tuple of surfaces and detector elements
  std::tuple<Surfaces, DetectorElements> convertGeant4Volumes() const;
};

}  // namespace ActsExamples::Geant4

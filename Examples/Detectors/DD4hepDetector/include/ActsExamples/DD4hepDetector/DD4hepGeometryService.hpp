// This file is part of the Acts project.
//
// Copyright (C) 2017 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "ActsExamples/Framework/BareService.hpp"
#include "ActsExamples/Framework/ProcessCode.hpp"
#include <Acts/Definitions/Units.hpp>
#include <Acts/Geometry/TrackingGeometry.hpp>
#include <Acts/Material/IMaterialDecorator.hpp>
#include <Acts/Utilities/BinningType.hpp>
#include <Acts/Utilities/Logger.hpp>

#include <functional>
#include <memory>

#include <DD4hep/DetElement.h>
#include <DD4hep/Detector.h>
#include <TGeoNode.h>

namespace ActsExamples {
namespace DD4hep {

void sortFCChhDetElements(std::vector<dd4hep::DetElement>& det);

/// @class DD4hepGeometryService
///
/// @brief service creating geometries from dd4hep input
///
/// The DD4hepGeometryService creates the DD4hep, the TGeo and the ACTS
/// TrackingGeometry
/// from DD4hep xml input. The geometries are created only on demand.
class DD4hepGeometryService final : public BareService {
 public:
  struct Config {
    /// Log level for the geometry service.
    Acts::Logging::Level logLevel = Acts::Logging::Level::INFO;
    /// Log level for DD4hep itself
    Acts::Logging::Level dd4hepLogLevel = Acts::Logging::Level::INFO;
    /// XML-file with the detector description
    std::vector<std::string> xmlFileNames;
    /// The name of the service
    std::string name;
    /// Binningtype in phi
    Acts::BinningType bTypePhi = Acts::equidistant;
    /// Binningtype in r
    Acts::BinningType bTypeR = Acts::arbitrary;
    /// Binningtype in z
    Acts::BinningType bTypeZ = Acts::equidistant;
    /// The tolerance added to the geometrical extension in r
    /// of the layers contained to build the volume envelope around
    /// @note this parameter only needs to be set if the volumes containing
    /// the
    /// layers (e.g. barrel, endcap volumes) have no specific shape
    /// (assemblies)
    double envelopeR = 1 * Acts::UnitConstants::mm;
    /// The tolerance added to the geometrical extension in z
    /// of the layers contained to build the volume envelope around
    /// @note this parameter only needs to be set if the volumes containing
    /// the layers (e.g. barrel, endcap volumes) have no specific shape
    /// (assemblies)
    double envelopeZ = 1 * Acts::UnitConstants::mm;
    double defaultLayerThickness = 10e-10;
    std::function<void(std::vector<dd4hep::DetElement>& detectors)>
        sortDetectors = sortFCChhDetElements;
    /// Material decorator
    std::shared_ptr<const Acts::IMaterialDecorator> matDecorator = nullptr;
  };

  DD4hepGeometryService(const Config& cfg);
  ~DD4hepGeometryService() final override;

  /// Interface method to access the DD4hep geometry
  /// @return The world DD4hep DetElement
  dd4hep::DetElement dd4hepGeometry();

  /// Interface method to Access the TGeo geometry
  /// @return The world TGeoNode (physical volume)
  TGeoNode* tgeoGeometry();

  /// Interface method to access to the interface of the DD4hep geometry
  dd4hep::Detector* lcdd();

  /// Interface method to access the ACTS TrackingGeometry
  ///
  /// @param gctx is the geometry context object
  std::unique_ptr<const Acts::TrackingGeometry> trackingGeometry(
      const Acts::GeometryContext& gctx);

 private:
  /// Private method to initiate building of the DD4hep geometry
  ActsExamples::ProcessCode buildDD4hepGeometry();

  /// Private method to initiate building of the ACTS tracking geometry
  ActsExamples::ProcessCode buildTrackingGeometry(
      const Acts::GeometryContext& gctx);

  /// The config class
  Config m_cfg;
  /// Pointer to the interface to the DD4hep geometry
  dd4hep::Detector* m_lcdd = nullptr;
  /// The world DD4hep DetElement
  dd4hep::DetElement m_dd4hepGeometry;
  /// The ACTS TrackingGeometry
  std::unique_ptr<const Acts::TrackingGeometry> m_trackingGeometry;
};

}  // namespace DD4hep
}  // namespace ActsExamples

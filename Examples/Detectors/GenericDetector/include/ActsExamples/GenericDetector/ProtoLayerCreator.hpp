// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Geometry/ApproachDescriptor.hpp"
#include "Acts/Geometry/DetectorElementBase.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/ProtoLayer.hpp"
#include "Acts/Material/ISurfaceMaterial.hpp"
#include "Acts/Surfaces/PlanarBounds.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/Surfaces/RectangleBounds.hpp"
#include "Acts/Surfaces/TrapezoidBounds.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/GenericDetector/GenericDetectorElement.hpp"

#include <iostream>

namespace Acts {
class LayerCreator;
class Surface;
class DetecorElementBase;
}  // namespace Acts

namespace ActsExamples::Generic {

struct ProtoLayerSurfaces {
  Acts::MutableProtoLayer protoLayer;
  std::vector<std::shared_ptr<Acts::Surface>> surfaces;
  std::size_t bins0;
  std::size_t bins1;
};

/// @class ProtoLayerCreatorT
///
/// The ProtoLayerCreatorT is the first setp in creating a geometry
/// from code input, it creates the ProtoLayer and returns the
/// created detector elements for the DetectorStore emulation
class ProtoLayerCreator {
 public:
  using DetectorElementFactory =
      std::function<std::shared_ptr<GenericDetectorElement>(
          std::shared_ptr<const Acts::Transform3>,
          std::shared_ptr<const Acts::PlanarBounds>, double,
          std::shared_ptr<const Acts::ISurfaceMaterial>)>;

  /// @struct Config
  ///
  /// Nested configuration struct for the ProtoLayerCreatorT
  struct Config {
    /// a single parameter for the approach surface envelope
    double approachSurfaceEnvelope = 0.5;
    /// central layer specification
    /// bin multipliers in rphi,z for finer module binning
    std::pair<int, int> centralLayerBinMultipliers;
    /// layer radii for the sensitive layers
    std::vector<double> centralLayerRadii;
    /// the (additional) layer envelope in R/Z
    std::vector<std::pair<double, double>> centralLayerEnvelopes;
    /// the binning schema: nPhi x nZ
    std::vector<std::pair<int, int>> centralModuleBinningSchema;
    /// the module center positions
    std::vector<std::vector<Acts::Vector3>> centralModulePositions;
    /// the module tilt for this layer
    std::vector<double> centralModuleTiltPhi;
    /// the module bounds: local x
    std::vector<double> centralModuleHalfX;
    /// the module bounds: local y
    std::vector<double> centralModuleHalfY;
    /// the module bounds: local z -> thickness
    std::vector<double> centralModuleThickness;
    /// the module material
    std::vector<std::shared_ptr<const Acts::ISurfaceMaterial>>
        centralModuleMaterial;
    /// the module front side stereo (if exists)
    std::vector<double> centralModuleFrontsideStereo;
    /// the module back side stereo (if exists)
    std::vector<double> centralModuleBacksideStereo;
    /// the module gap between frontside and backside
    std::vector<double> centralModuleBacksideGap;

    /// the layers at p/e side
    /// bin multipliers in r,phi for finer module binning
    std::pair<int, int> posnegLayerBinMultipliers;
    /// layer positions in Z
    std::vector<double> posnegLayerPositionsZ;
    /// the envelope definitions
    std::vector<double> posnegLayerEnvelopeR;
    /// the module center positions
    std::vector<std::vector<std::vector<Acts::Vector3>>> posnegModulePositions;
    /// the phi binning
    std::vector<std::vector<std::size_t>> posnegModulePhiBins;
    /// the module bounds: min halfx
    std::vector<std::vector<double>> posnegModuleMinHalfX;
    /// the module bounds: max halfx
    std::vector<std::vector<double>> posnegModuleMaxHalfX;
    /// the module bounds: local y
    std::vector<std::vector<double>> posnegModuleHalfY;
    /// the module bounds: local z -> thickness
    std::vector<std::vector<double>> posnegModuleThickness;
    /// the module material
    std::vector<std::vector<std::shared_ptr<const Acts::ISurfaceMaterial>>>
        posnegModuleMaterial;
    /// the module front side stereo (if exists)
    std::vector<std::vector<double>> posnegModuleFrontsideStereo;
    /// the module back side stereo (if exists)
    std::vector<std::vector<double>> posnegModuleBacksideStereo;
    /// the module gap between frontside and backside
    std::vector<std::vector<double>> posnegModuleBacksideGap;

    /// The factory for the detector elements. This function should make sure
    /// that the lifetime of the detector element is managed by the creator.
    /// It's also responsible for creating the identifiers for the detector
    /// elements.
    DetectorElementFactory detectorElementFactory;
  };

  /// Constructor
  /// @param cfg is the configuration class
  /// @param logger is the logging class for screen output
  explicit ProtoLayerCreator(const Config& cfg,
                             std::unique_ptr<const Acts::Logger> logger =
                                 Acts::getDefaultLogger("ProtoLayerCreator",
                                                        Acts::Logging::INFO));

  /// @brief construct the negative side layers
  /// @param gctx The geometry context for this construction call
  /// @param detectorStore The reference store for the detector elements
  /// @return the protolayers and surfaces on the negative detector side
  std::vector<ProtoLayerSurfaces> negativeProtoLayers(
      const Acts::GeometryContext& gctx) const;

  /// @brief construct the central layers
  /// @param gctx The geometry context for this construction call
  /// @param detectorStore The reference store for the detector elements
  /// @return the protolayers and surfaces on the central detector side
  std::vector<ProtoLayerSurfaces> centralProtoLayers(
      const Acts::GeometryContext& gctx) const;

  /// @brief construct the positive side layers
  /// @param gctx The geometry context for this construction call
  /// @param detectorStore The reference store for the detector elements
  /// @return the protolayers and surfaces on the  positive detector side
  std::vector<ProtoLayerSurfaces> positiveProtoLayers(
      const Acts::GeometryContext& gctx) const;

  /// @brief get the configuration
  /// @return the configuration
  const Config& config() const { return m_cfg; }

 private:
  /// @brief private helper method to create the proto layers on the
  /// left respectively right side
  /// @param gctx The geometry context for this construction call
  /// @param detectorStore The reference store for the detector elements
  /// @param side is the indiciator whether to build on negative/positive
  /// @return the protolayers and surfaces on the neg/pos detector side
  std::vector<ProtoLayerSurfaces> createProtoLayers(
      const Acts::GeometryContext& gctx, int side) const;

  /// Configuration member
  Config m_cfg;

  /// the logging instance
  std::unique_ptr<const Acts::Logger> m_logger;

  /// Private access to the logging instance
  const Acts::Logger& logger() const { return *m_logger; }
};

}  // namespace ActsExamples::Generic

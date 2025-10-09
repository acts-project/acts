// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/Units.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/LayerArrayCreator.hpp"
#include "Acts/Geometry/LayerCreator.hpp"
#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/Geometry/TrackingVolume.hpp"
#include "Acts/Tests/CommonHelpers/DetectorElementStub.hpp"
#include "Acts/Tests/CommonHelpers/PredefinedMaterials.hpp"
#include "Acts/Utilities/Logger.hpp"

#include <vector>

namespace Acts::Test {

struct CylindricalTrackingGeometry {
  constexpr static double kBeamPipeRadius = 19. * Acts::UnitConstants::mm;
  constexpr static double kBeamPipeHalfLengthZ =
      1000. * Acts::UnitConstants::mm;
  constexpr static double kBeamPipeThickness = 0.8 * Acts::UnitConstants::mm;

  inline static const std::vector<double> kLayerRadii = {32., 72., 116., 172.};
  inline static const std::vector<std::pair<int, int>> kLayerBinning = {
      {16, 14}, {32, 14}, {52, 14}, {78, 14}};
  inline static const std::vector<double> kModuleTiltPhi = {0.145, 0.145, 0.145,
                                                            0.145};
  inline static const std::vector<double> kModuleHalfX = {8.4, 8.4, 8.4, 8.4};
  inline static const std::vector<double> kModuleHalfY = {36., 36., 36., 36.};
  inline static const std::vector<double> kModuleThickness = {0.15, 0.15, 0.15,
                                                              0.15};

  std::reference_wrapper<const GeometryContext> geoContext;
  bool gen3 = false;

  /// Only allowed constructor with reference wrapper
  explicit CylindricalTrackingGeometry(const GeometryContext& gctx,
                                       bool gen3_ = false)
      : geoContext(gctx), gen3(gen3_) {}

  using DetectorStore = std::vector<std::unique_ptr<const DetectorElementStub>>;

  /// The detector store for memory management
  DetectorStore detectorStore = {};

  /// Generator of surfaces for a ring
  ///
  /// @param detStore The DetectorStore for storing the modules
  /// @param moduleHalfXminY The half length in X (at Y min) of the module
  /// @param moduleHalfXmaxY The half length in X (at Y max) of the module
  /// @param moduleHalfY The half length in Y of the module
  /// @param moduleThickness The module thickness
  /// @param moduleTilt The tilt out of the plane for discs
  /// @param ringRadius The central radius of the ring
  /// @param ringZ The z position of the ring
  /// @param zStagger The z offset of phi modules
  /// @param nPhi The number of phi modules
  ///
  /// @return A vector of Surfaces
  std::vector<Surface*> surfacesRing(DetectorStore& detStore,
                                     double moduleHalfXminY,
                                     double moduleHalfXmaxY, double moduleHalfY,
                                     double moduleThickness, double moduleTilt,
                                     double ringRadius, double ringZ,
                                     double zStagger, int nPhi);

  /// Generator of surfaces for a cylindrical layer
  ///
  /// @param detStore The DetectorStore for storing the modules
  /// @param moduleHalfX The half length in X of the module
  /// @param moduleHalfY The half length in Y of the module
  /// @param moduleThickness The module thickness
  /// @param moduleTilePhi The tilt in phi direction of the module
  /// @param layerRadius The radius of the cylindrical layer
  /// @param radialStagger The radial delta of modules next in z
  /// @param longitudinalOverlap The z overlap of modules next in z
  /// @param binningSchema The number of bins in phi/z
  ///
  /// @return A vector of Surfaces
  std::vector<Surface*> surfacesCylinder(
      DetectorStore& detStore, double moduleHalfX, double moduleHalfY,
      double moduleThickness, double moduleTiltPhi, double layerRadius,
      double radialStagger, double longitudinalOverlap,
      const std::pair<int, int>& binningSchema);

  /// Helper method for cylinder layer
  /// create the positions for module surfaces on a cylinder
  std::vector<Vector3> modulePositionsCylinder(
      double radius, double zStagger, double moduleHalfLength, double lOverlap,
      const std::pair<int, int>& binningSchema);

  // @brief Call operator for the creation method of the tracking geometry
  std::shared_ptr<TrackingGeometry> operator()(
      const Logger& logger = getDummyLogger());

 private:
  std::shared_ptr<TrackingGeometry> buildGen1(const Logger& logger);
  std::shared_ptr<TrackingGeometry> buildGen3(const Logger& logger);
};

}  // namespace Acts::Test

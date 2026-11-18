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
#include "Acts/Utilities/Logger.hpp"
#include "ActsTests/CommonHelpers/DetectorElementStub.hpp"

#include <vector>

namespace ActsTests {

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

  /// Gen3 endcaps (gen3Endcaps): |z| of the discs, mirrored at z = 0
  inline static const std::vector<double> kDiscZ = {560., 660., 780., 920.};
  /// A ring of trapezoid modules of a disc, local y radial
  struct DiscRing {
    double radius;
    double halfXminY;
    double halfXmaxY;
    double halfY;
    int nPhi;
    /// z offset from the disc position, away from z = 0
    double zOffset;
  };
  /// Inner ring r 29-81 mm, outer ring r 80-176 mm
  inline static const std::vector<DiscRing> kDiscRings = {
      {55., 4.2, 11.5, 26., 24, -2.}, {128., 6.8, 14.5, 48., 40, 2.}};
  /// z offset of every second module in phi, away from z = 0
  constexpr static double kDiscModuleStagger = 1. * Acts::UnitConstants::mm;

  std::reference_wrapper<const Acts::GeometryContext> geoContext;
  bool gen3 = false;
  /// Gen3 only: fill the radial gaps between layers with gap volumes instead
  /// of extending the layer volumes to the midpoint (default: midpoint)
  bool gen3GapVolumes = false;
  /// Gen3 only: radial envelope of the layer volumes around the modules
  double gen3LayerEnvelopeR = 5. * Acts::UnitConstants::mm;
  /// Gen3 only: give the layer volumes their portals through
  /// CylinderNavigationPolicy (as ODD, GenericDetector and ATLAS ITk do)
  /// instead of TryAllNavigationPolicy{.sensitives = false}
  bool gen3CylinderPortalPolicy = false;
  /// Gen3 only: add two endcaps of four disc layers (kDiscZ, kDiscRings) next
  /// to the barrel. The barrel layers then sit in their own container; the
  /// disc layers take gen3LayerEnvelopeR as their z envelope.
  bool gen3Endcaps = false;
  /// Gen3 endcaps only: stretch the disc layer volumes until they meet
  /// (Midpoint attachment, Expand resize) instead of gap volumes between them
  bool gen3ThickEndcaps = false;
  /// Gen3 only: passive material on the beampipe and on the layer portals
  /// (1.5 mm Si each). False leaves the module material (kModuleThickness of
  /// Si) as the only material.
  bool gen3PassiveMaterial = true;

  /// Only allowed constructor with reference wrapper
  explicit CylindricalTrackingGeometry(const Acts::GeometryContext& gctx,
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
  std::vector<Acts::Surface*> surfacesRing(
      DetectorStore& detStore, double moduleHalfXminY, double moduleHalfXmaxY,
      double moduleHalfY, double moduleThickness, double moduleTilt,
      double ringRadius, double ringZ, double zStagger, int nPhi);

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
  std::vector<Acts::Surface*> surfacesCylinder(
      DetectorStore& detStore, double moduleHalfX, double moduleHalfY,
      double moduleThickness, double moduleTiltPhi, double layerRadius,
      double radialStagger, double longitudinalOverlap,
      const std::pair<int, int>& binningSchema);

  /// Helper method for cylinder layer
  /// create the positions for module surfaces on a cylinder
  std::vector<Acts::Vector3> modulePositionsCylinder(
      double radius, double zStagger, double moduleHalfLength, double lOverlap,
      const std::pair<int, int>& binningSchema);

  // @brief Call operator for the creation method of the tracking geometry
  std::shared_ptr<Acts::TrackingGeometry> operator()(
      const Acts::Logger& logger = Acts::getDummyLogger());

 private:
  std::shared_ptr<Acts::TrackingGeometry> buildGen1(const Acts::Logger& logger);
  std::shared_ptr<Acts::TrackingGeometry> buildGen3(const Acts::Logger& logger);
};

}  // namespace ActsTests

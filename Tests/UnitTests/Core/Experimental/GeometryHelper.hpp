// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/Units.hpp"
#include "Acts/Experimental/CylindricalContainerHelper.hpp"
#include "Acts/Experimental/DetectorVolume.hpp"
#include "Acts/Experimental/SurfaceLinks.hpp"
#include "Acts/Geometry/CylinderVolumeBounds.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/Surfaces/RectangleBounds.hpp"
#include "Acts/Surfaces/TrapezoidBounds.hpp"
#include "Acts/Tests/CommonHelpers/CylindricalTrackingGeometry.hpp"

namespace Acts {
namespace Test {

/// Generator of surfaces for a cylindrical layer
///
/// @param moduleHalfX The half lenght in X of the module
/// @param moduleHalfY The half lenght in Y of the module
/// @param moduleThickness The module thickness
/// @param moduleTilePhi The tilt in phi direction of the module
/// @param layerRadius The radius fo the cylindrical layer
/// @param radialStagger The radial delta of modules next in z
/// @param longitudinalOverlap The z overlap of modules next in z
/// @param binningSchema The number of bins in phi/z
///
/// @return A vector of Surfaces
std::vector<std::shared_ptr<Surface>> surfacesCylinder(
    ActsScalar moduleHalfX, ActsScalar moduleHalfY, ActsScalar moduleTiltPhi,
    ActsScalar layerRadius, ActsScalar radialStagger,
    ActsScalar longitudinalOverlap, const std::pair<int, int>& binningSchema) {
  // Prepare the return vector
  std::vector<std::shared_ptr<Surface>> layerSurfaces;

  // The rectangle bounds for all modules
  auto mBounds = std::make_shared<RectangleBounds>(moduleHalfX, moduleHalfY);

  // Create the module centers
  auto moduleCenters = CylindricalTrackingGeometry::modulePositionsCylinder(
      layerRadius, radialStagger, moduleHalfY, longitudinalOverlap,
      binningSchema);

  for (auto& mCenter : moduleCenters) {
    // The association transform
    ActsScalar modulePhi = VectorHelpers::phi(mCenter);
    // Local z axis is the normal vector
    Vector3 moduleLocalZ(cos(modulePhi + moduleTiltPhi),
                         sin(modulePhi + moduleTiltPhi), 0.);
    // Local y axis is the global z axis
    Vector3 moduleLocalY(0., 0., 1);
    // Local x axis the normal to local y,z
    Vector3 moduleLocalX(-sin(modulePhi + moduleTiltPhi),
                         cos(modulePhi + moduleTiltPhi), 0.);
    // Create the RotationMatrix
    RotationMatrix3 moduleRotation;
    moduleRotation.col(0) = moduleLocalX;
    moduleRotation.col(1) = moduleLocalY;
    moduleRotation.col(2) = moduleLocalZ;
    // Get the moduleTransform
    auto mModuleTransform = Transform3(Translation3(mCenter) * moduleRotation);

    layerSurfaces.push_back(
        Surface::makeShared<PlaneSurface>(mModuleTransform, mBounds));
  }
  return layerSurfaces;
}

/// Generator of surfaces for a ring
///
/// @param moduleHalfXminY The half lenght in X (at Y min) of the module
/// @param moduleHalfXmaxY The half lenght in X (at Y max) of the module
/// @param moduleHalfY The half lenght in Y of the module
/// @param moduleTilt The tilt out of the plane for discs
/// @param ringRadius The central radius of the ring
/// @param ringZ The z position of the ring
/// @param zStagger The z offset of phi moudles
/// @param nPhi The number of phi modules
///
/// @return A vector of Surfaces
std::vector<std::shared_ptr<Surface>> surfacesRing(
    ActsScalar moduleHalfXminY, ActsScalar moudleHalfXmaxY,
    ActsScalar moduleHalfY, ActsScalar moduleTilt, ActsScalar ringRadius,
    ActsScalar ringZ, ActsScalar zStagger, int nPhi) {
  // Prepare the return vector
  std::vector<std::shared_ptr<Surface>> layerSurfaces;

  // The rectangle/trapezoid bounds for all modules
  std::shared_ptr<PlanarBounds> mBounds = nullptr;
  if (moduleHalfXminY == moudleHalfXmaxY) {
    mBounds = std::make_shared<RectangleBounds>(moduleHalfXminY, moduleHalfY);
  } else {
    mBounds = std::make_shared<TrapezoidBounds>(moduleHalfXminY,
                                                moudleHalfXmaxY, moduleHalfY);
  }

  ActsScalar phiStep = 2 * M_PI / nPhi;

  for (int im = 0; im < nPhi; ++im) {
    // Get the moduleTransform
    ActsScalar phi = -M_PI + im * phiStep;
    auto mModuleTransform = Transform3(
        Translation3(ringRadius * std::cos(phi), ringRadius * std::sin(phi),
                     ringZ + (im % 2) * zStagger) *
        AngleAxis3(phi - 0.5 * M_PI, Vector3::UnitZ()) *
        AngleAxis3(moduleTilt, Vector3::UnitY()));

    // Create the detector element
    auto detSurface =
        Surface::makeShared<PlaneSurface>(mModuleTransform, mBounds);
    layerSurfaces.push_back(detSurface);
  }

  return layerSurfaces;
}

/// Create a single layer Volume -
///
/// @param volumeMinR minimal radius of volume
/// @param volumeMmaxR maximal radius of volume
/// @param volumeHalfZ half length in z of volume
/// @param moduleHalfX The half lenght in X of the module
/// @param moduleHalfY The half lenght in Y of the module
/// @param moduleThickness The module thickness
/// @param moduleTilePhi The tilt in phi direction of the module
/// @param layerRadius The radius fo the cylindrical layer
/// @param radialStagger The radial delta of modules next in z
/// @param longitudinalOverlap The z overlap of modules next in z
/// @param binningSchema The number of bins in phi/z
///
/// @return a container volume for a barrel
std::shared_ptr<DetectorVolume> createBarrelVolume(
    ActsScalar volumeMinR = 25., ActsScalar volumeMaxR = 40.,
    ActsScalar volumeHalfZ = 500.,
    const std::string& volumeName = "SingleLayerVolume",
    ActsScalar moduleHalfX = 8.4, ActsScalar moduleHalfY = 36.,
    ActsScalar moduleTiltPhi = 0.145, ActsScalar layerRadius = 32.,
    ActsScalar radialStagger = 2., ActsScalar longitudinalOverlap = 5.,
    const std::pair<int, int>& binningSchema = {16, 14}) {
  // Generate the volume surfaces
  std::vector<std::shared_ptr<Surface>> volumeSurfaces =
      surfacesCylinder(moduleHalfX, moduleHalfY, moduleTiltPhi, layerRadius,
                       radialStagger, longitudinalOverlap, binningSchema);

  // Create the volume bounds
  auto volumeBounds = std::make_unique<CylinderVolumeBounds>(
      volumeMinR, volumeMaxR, volumeHalfZ);

  SurfaceLinks volumeSurfaceLinks = AllSurfaces{};

  std::vector<SurfaceLinks> portalSurfaceLinks = {AllSurfaces{}, AllSurfaces{},
                                                  AllSurfaces{}, AllSurfaces{}};

  return DetectorVolume::makeShared(
      Transform3::Identity(), std::move(volumeBounds),
      std::move(volumeSurfaces), std::move(volumeSurfaceLinks),
      std::move(portalSurfaceLinks), volumeName);
}

/// @param volumeMinR The volume rmin
/// @param volumeMaxR The volume rmax
/// @param volumeMinZ The volume zmin
/// @param volumeMaxZ The volume zmax
/// @param moduleHalfXminY The half lenght in X (at Y min) of the module
/// @param moduleHalfXmaxY The half lenght in X (at Y max) of the module
/// @param moduleHalfY The half lenght in Y of the module
/// @param moduleTilt The tilt out of the plane for discs
/// @param ringRadius The central radius of the ring
/// @param ringZ The z position of the ring
/// @param zStagger The z offset of phi moudles
/// @param nPhi The number of phi modules
///
/// @return an endcap volume
std::shared_ptr<DetectorVolume> createEndcapVolume(
    ActsScalar volumeMinR, ActsScalar volumeMaxR, ActsScalar volumeMinZ,
    ActsScalar volumeMaxZ, ActsScalar moduleHalfXminY = 8.4,
    ActsScalar moudleHalfXmaxY = 12.4, ActsScalar moduleHalfY = 32.,
    ActsScalar moduleTilt = 0., ActsScalar ringRadius = 40.,
    ActsScalar ringZ = 600., ActsScalar zStagger = 2, int nPhi = 40,
    const std::string& volumeName = "Unnamed") {
  auto volumeSurfaces =
      surfacesRing(moduleHalfXminY, moudleHalfXmaxY, moduleHalfY, moduleTilt,
                   ringRadius, ringZ, zStagger, nPhi);

  // Create the volume bounds
  auto volumeBounds = std::make_unique<CylinderVolumeBounds>(
      volumeMinR, volumeMaxR, 0.5 * (volumeMaxZ - volumeMinZ));

  SurfaceLinks volumeSurfaceLinks = AllSurfaces{};

  std::vector<SurfaceLinks> portalSurfaceLinks = {AllSurfaces{}, AllSurfaces{},
                                                  AllSurfaces{}};

  if (volumeMinR > 0.) {
    portalSurfaceLinks.push_back(AllSurfaces{});
  }

  auto volumeTransform = Transform3::Identity();
  volumeTransform.pretranslate(Vector3(0., 0., 3.));

  return DetectorVolume::makeShared(
      volumeTransform, std::move(volumeBounds), std::move(volumeSurfaces),
      std::move(volumeSurfaceLinks), std::move(portalSurfaceLinks), volumeName);
}

/// Helper method to create a central detector
///
/// Keep track of the central half length with
/// @param detectorHalfZ
///
/// @return a central detector volume
std::shared_ptr<DetectorVolume> createCentralDetector(
    ActsScalar detectorHalfZ = 500.) {
  // Create the volume bounds
  ActsScalar beamPipeR = 27.;
  // Beam pipe volume
  auto beamPipeBounds =
      std::make_unique<CylinderVolumeBounds>(0., beamPipeR, detectorHalfZ);
  auto beamPipe = DetectorVolume::makeShared(
      Transform3::Identity(), std::move(beamPipeBounds), "BeamPipe");
  // First layer
  ActsScalar firstLayerOuterR = 38.;
  auto firstLayer = createBarrelVolume(beamPipeR, firstLayerOuterR,
                                       detectorHalfZ, "BarrelLayer0");
  // First gap
  ActsScalar secondLayerInnerR = 64.;
  auto firstGapBounds = std::make_unique<CylinderVolumeBounds>(
      firstLayerOuterR, secondLayerInnerR, detectorHalfZ);
  auto firstGap = DetectorVolume::makeShared(
      Transform3::Identity(), std::move(firstGapBounds), "BarrelGap0");
  // Second layer
  ActsScalar secondLayerOuterR = 80.;
  auto secondLayer = createBarrelVolume(secondLayerInnerR, secondLayerOuterR,
                                        detectorHalfZ, "BarrelLayer1", 8.4, 36.,
                                        0.145, 72., 2., 5., {32, 14});

  // The volumes in R
  std::vector<std::shared_ptr<DetectorVolume>> barrelVolumes = {
      beamPipe, firstLayer, firstGap, secondLayer};

  // Return the container in R
  return CylindricalContainerHelper::containerInR(std::move(barrelVolumes),
                                                  "BarrelWithTwoLayers");
}

}  // namespace Test
}  // namespace Acts
// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Geometry/AbstractVolume.hpp"
#include "Acts/Geometry/ConeVolumeBounds.hpp"
#include "Acts/Geometry/CuboidVolumeBounds.hpp"
#include "Acts/Geometry/CylinderVolumeBounds.hpp"
#include "Acts/Geometry/GenericCuboidVolumeBounds.hpp"
#include "Acts/Geometry/TrapezoidVolumeBounds.hpp"
#include "Acts/Visualization/GeometryView3D.hpp"
#include "Acts/Visualization/IVisualization3D.hpp"

#include <fstream>
#include <sstream>
#include <string>

namespace Acts {

namespace VolumeView3DTest {

/// Helper method to visualize all types of surfaces
///
/// @param helper The visualization helper
/// @param triangulate The directive whether to create triangular meshes
/// @param tag The test tag (mode) identification
///
/// @return a string containing all written characters

static inline std::string run(IVisualization3D& helper, bool triangulate,
                              const std::string& tag) {
  auto gctx = GeometryContext();
  auto identity = Transform3::Identity();
  std::stringstream cStream;

  double halfPhiSector = M_PI / 4.;

  ViewConfig vConfig = s_viewVolume;
  vConfig.triangulate = triangulate;

  // ---------------------------------------------------
  // Cuboid surface section
  auto box = std::make_shared<CuboidVolumeBounds>(4., 3., 6.);
  auto cuboid = std::make_shared<AbstractVolume>(identity, box);
  GeometryView3D::drawVolume(helper, *cuboid, gctx, Transform3::Identity(),
                             vConfig);
  helper.write(std::string("Volumes_CuboidVolume") + tag);
  helper.write(cStream);
  helper.clear();

  //----------------------------------------------------
  // Cone volume section
  // Single solid Cone
  auto solidCone =
      std::make_shared<ConeVolumeBounds>(0., 0., 0.45, 5., 5., 0., M_PI);
  auto cone = std::make_shared<AbstractVolume>(identity, solidCone);
  GeometryView3D::drawVolume(helper, *cone, gctx, Transform3::Identity(),
                             vConfig);
  helper.write("Volumes_ConeVolumeSolid");
  helper.write(cStream);
  helper.clear();

  // Single solid Cone - with cut off
  auto cutOffCone =
      std::make_shared<ConeVolumeBounds>(0., 0., 0.45, 8., 5., 0., M_PI);
  cone = std::make_shared<AbstractVolume>(identity, cutOffCone);
  GeometryView3D::drawVolume(helper, *cone, gctx, Transform3::Identity(),
                             vConfig);
  helper.write("Volumes_ConeVolumeSolidCutOff");
  helper.write(cStream);
  helper.clear();

  // Cone - Cone inlay
  auto cutOffHollowCone =
      std::make_shared<ConeVolumeBounds>(0.35, 7., 0.45, 8., 5, 0., M_PI);
  cone = std::make_shared<AbstractVolume>(identity, cutOffHollowCone);
  GeometryView3D::drawVolume(helper, *cone, gctx, Transform3::Identity(),
                             vConfig);
  helper.write("Volumes_ConeVolumeConeCone");
  helper.write(cStream);
  helper.clear();

  // Sectoral Cone - Cone inlay
  auto cutOffHollowSectoralCone =
      std::make_shared<ConeVolumeBounds>(0.35, 7., 0.45, 8., 5., 0., 0.456);
  cone = std::make_shared<AbstractVolume>(identity, cutOffHollowSectoralCone);
  GeometryView3D::drawVolume(helper, *cone, gctx, Transform3::Identity(),
                             vConfig);
  helper.write("Volumes_ConeVolumeConeConeSectoral");
  helper.write(cStream);
  helper.clear();

  // Single Hollow Cone - cylindrical inlay
  auto cutOffHollowCylCone =
      std::make_shared<ConeVolumeBounds>(1., 0.45, 8., 5., 0., M_PI);
  cone = std::make_shared<AbstractVolume>(identity, cutOffHollowCylCone);
  GeometryView3D::drawVolume(helper, *cone, gctx, Transform3::Identity(),
                             vConfig);
  helper.write("Volumes_ConeVolumeConeCylinder");
  helper.write(cStream);
  helper.clear();

  // Single Hollow Cylinder - Cone inlay
  auto cutOffHollowConeCyl =
      std::make_shared<ConeVolumeBounds>(12., 0.35, 7., 5., 0., M_PI);
  cone = std::make_shared<AbstractVolume>(identity, cutOffHollowConeCyl);
  GeometryView3D::drawVolume(helper, *cone, gctx, Transform3::Identity(),
                             vConfig);
  helper.write("Volumes_ConeVolumeCylinderCone");
  helper.write(cStream);
  helper.clear();

  //----------------------------------------------------
  // Cylinder volume section
  double cylinderInnerR = 1.;
  double cylinderOuterR = 5.;
  double cylinderHalfZ = 10.;

  auto fullCylinder =
      std::make_shared<CylinderVolumeBounds>(0., cylinderOuterR, cylinderHalfZ);
  auto cylinder = std::make_shared<AbstractVolume>(identity, fullCylinder);
  GeometryView3D::drawVolume(helper, *cylinder, gctx, Transform3::Identity(),
                             vConfig);
  helper.write("Volumes_CylinderVolumeFull");
  helper.write(cStream);
  helper.clear();

  auto tubeCylinder = std::make_shared<CylinderVolumeBounds>(
      cylinderInnerR, cylinderOuterR, cylinderHalfZ);
  cylinder = std::make_shared<AbstractVolume>(identity, tubeCylinder);
  GeometryView3D::drawVolume(helper, *cylinder, gctx, Transform3::Identity(),
                             vConfig);
  helper.write("Volumes_CylinderVolumeTube");
  helper.write(cStream);
  helper.clear();

  tubeCylinder = std::make_shared<CylinderVolumeBounds>(
      cylinderInnerR, cylinderOuterR, cylinderHalfZ, halfPhiSector);
  cylinder = std::make_shared<AbstractVolume>(identity, tubeCylinder);
  GeometryView3D::drawVolume(helper, *cylinder, gctx, Transform3::Identity(),
                             vConfig);
  helper.write("Volumes_CylinderVolumeTubeSector");
  helper.write(cStream);
  helper.clear();

  //----------------------------------------------------
  // Trapezoid volume section
  std::array<Vector3, 8> vertices;
  vertices = {{{0, 0, 0},
               {2, 0, 0},
               {2, 1, 0},
               {0, 1, 0},
               {0, 0, 1},
               {2, 0, 1},
               {2, 1, 1},
               {0, 1, 1}}};
  auto genericCuboid = std::make_shared<GenericCuboidVolumeBounds>(vertices);
  auto generic = std::make_shared<AbstractVolume>(identity, genericCuboid);
  GeometryView3D::drawVolume(helper, *generic, gctx, Transform3::Identity(),
                             vConfig);
  helper.write("Volumes_GenericCuboidVolume");
  helper.write(cStream);
  helper.clear();

  //----------------------------------------------------
  // Trapezoid volume section
  auto trapezoid = std::make_shared<TrapezoidVolumeBounds>(2., 4., 5., 6.);
  auto trapezoidVolume = std::make_shared<AbstractVolume>(identity, trapezoid);
  GeometryView3D::drawVolume(helper, *trapezoidVolume, gctx,
                             Transform3::Identity(), vConfig);
  helper.write("Volumes_TrapezoidVolume");
  helper.write(cStream);
  helper.clear();

  return cStream.str();
}

}  // namespace VolumeView3DTest
}  // namespace Acts

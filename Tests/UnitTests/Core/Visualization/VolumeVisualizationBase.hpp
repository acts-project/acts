// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Visualization/GeometryVisualization.hpp"
#include "Acts/Visualization/IVisualization.hpp"

#include "Acts/Geometry/AbstractVolume.hpp"
#include "Acts/Geometry/ConeVolumeBounds.hpp"
#include "Acts/Geometry/CuboidVolumeBounds.hpp"
#include "Acts/Geometry/CylinderVolumeBounds.hpp"
#include "Acts/Geometry/GenericCuboidVolumeBounds.hpp"
#include "Acts/Geometry/TrapezoidVolumeBounds.hpp"

#include <fstream>
#include <sstream>
#include <string>

namespace Acts {
namespace VolumeVisualization {
/// Helper method to visualiza all types of surfaces
///
/// @param helper The visualziation helper
/// @param triangulate The directive whether to create triangular meshes
/// @param tag The test tag (mode) identification
///
/// @return a string containing all written caracters

static inline std::string test(IVisualization& helper, bool triangulate,
                               const std::string& tag) {
  auto gctx = GeometryContext();
  auto identity = std::make_shared<Transform3D>(Transform3D::Identity());
  std::ofstream stream;
  std::stringstream cStream;

  double halfPhiSector = M_PI / 4.;

  /// Helper method to prepare the streams & helpers
  /// @param path is the file path
  /// @param clear ist he indicator to clear the helper
  auto write = [&](const std::string& path, bool clear = true) -> void {
    std::string wpath = path + tag;
    helper.write(wpath);
    helper.write(cStream);
    if (clear) {
      helper.clear();
    }
  };

  // ---------------------------------------------------
  // Cuboid surface section
  IVisualization::ColorType boxColor = {0, 0, 255};

  auto box = std::make_shared<CuboidVolumeBounds>(4., 3., 6.);
  auto cuboid = std::make_shared<AbstractVolume>(identity, box);
  GeometryVisualization::drawVolume(helper, *cuboid, gctx,
                                    Transform3D::Identity(), 72, triangulate,
                                    boxColor);
  write("Volumes_CuboidVolume");

  //----------------------------------------------------
  // Cone volume section
  IVisualization::ColorType coneColor = {255, 171, 3};

  // Single solid Cone
  auto solidCone =
      std::make_shared<ConeVolumeBounds>(0., 0., 0.45, 5., 5., 0., M_PI);
  auto cone = std::make_shared<AbstractVolume>(identity, solidCone);
  GeometryVisualization::drawVolume(
      helper, *cone, gctx, Transform3D::Identity(), 72, triangulate, coneColor);
  write("Volumes_ConeVolumeSolid");

  // Single solid Cone - with cut off
  auto cutOffCone =
      std::make_shared<ConeVolumeBounds>(0., 0., 0.45, 8., 5., 0., M_PI);
  cone = std::make_shared<AbstractVolume>(identity, cutOffCone);
  GeometryVisualization::drawVolume(
      helper, *cone, gctx, Transform3D::Identity(), 72, triangulate, coneColor);
  write("Volumes_ConeVolumeSolidCutOff");

  // Cone - Cone inlay
  auto cutOffHollowCone =
      std::make_shared<ConeVolumeBounds>(0.35, 7., 0.45, 8., 5, 0., M_PI);
  cone = std::make_shared<AbstractVolume>(identity, cutOffHollowCone);
  GeometryVisualization::drawVolume(
      helper, *cone, gctx, Transform3D::Identity(), 72, triangulate, coneColor);
  write("Volumes_ConeVolumeConeCone");

  // Sectoral Cone - Cone inlay
  auto cutOffHollowSectoralCone =
      std::make_shared<ConeVolumeBounds>(0.35, 7., 0.45, 8., 5., 0., 0.456);
  cone = std::make_shared<AbstractVolume>(identity, cutOffHollowSectoralCone);
  GeometryVisualization::drawVolume(
      helper, *cone, gctx, Transform3D::Identity(), 72, triangulate, coneColor);
  write("Volumes_ConeVolumeConeConeSectoral");

  // Single Hollow Cone - cylindrical inlay
  auto cutOffHollowCylCone =
      std::make_shared<ConeVolumeBounds>(1., 0.45, 8., 5., 0., M_PI);
  cone = std::make_shared<AbstractVolume>(identity, cutOffHollowCylCone);
  GeometryVisualization::drawVolume(
      helper, *cone, gctx, Transform3D::Identity(), 72, triangulate, coneColor);
  write("Volumes_ConeVolumeConeCylinder");

  // Single Hollow Cylinder - Cone inlay
  auto cutOffHollowConeCyl =
      std::make_shared<ConeVolumeBounds>(12., 0.35, 7., 5., 0., M_PI);
  cone = std::make_shared<AbstractVolume>(identity, cutOffHollowConeCyl);
  GeometryVisualization::drawVolume(
      helper, *cone, gctx, Transform3D::Identity(), 72, triangulate, coneColor);
  write("Volumes_ConeVolumeCylinderCone");

  //----------------------------------------------------
  // Cylinder volume section
  IVisualization::ColorType cylinderColor = {0, 196, 252};

  double cylinderInnerR = 1.;
  double cylinderOuterR = 5.;
  double cylinderHalfZ = 10.;

  auto fullCylinder =
      std::make_shared<CylinderVolumeBounds>(0., cylinderOuterR, cylinderHalfZ);
  auto cylinder = std::make_shared<AbstractVolume>(identity, fullCylinder);
  GeometryVisualization::drawVolume(helper, *cylinder, gctx,
                                    Transform3D::Identity(), 72, triangulate,
                                    cylinderColor);
  write("Volumes_CylinderVolumeFull");

  auto tubeCylinder = std::make_shared<CylinderVolumeBounds>(
      cylinderInnerR, cylinderOuterR, cylinderHalfZ);
  cylinder = std::make_shared<AbstractVolume>(identity, tubeCylinder);
  GeometryVisualization::drawVolume(helper, *cylinder, gctx,
                                    Transform3D::Identity(), 72, triangulate,
                                    cylinderColor);
  write("Volumes_CylinderVolumeTube");

  tubeCylinder = std::make_shared<CylinderVolumeBounds>(
      cylinderInnerR, cylinderOuterR, cylinderHalfZ, halfPhiSector);
  cylinder = std::make_shared<AbstractVolume>(identity, tubeCylinder);
  GeometryVisualization::drawVolume(helper, *cylinder, gctx,
                                    Transform3D::Identity(), 72, triangulate,
                                    cylinderColor);
  write("Volumes_CylinderVolumeTubeSector");

  //----------------------------------------------------
  // Trapezoid volume section
  IVisualization::ColorType genericColor = {23, 196, 22};
  std::array<Vector3D, 8> vertices;
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
  GeometryVisualization::drawVolume(helper, *generic, gctx,
                                    Transform3D::Identity(), 72, triangulate,
                                    genericColor);
  write("Volumes_GenericCuboidVolume");

  //----------------------------------------------------
  // Trapezoid volume section
  IVisualization::ColorType trapezoidColor = {23, 196, 22};
  auto trapezoid = std::make_shared<TrapezoidVolumeBounds>(2., 4., 5., 6.);
  auto trapezoidVolume = std::make_shared<AbstractVolume>(identity, trapezoid);
  GeometryVisualization::drawVolume(helper, *trapezoidVolume, gctx,
                                    Transform3D::Identity(), 72, triangulate,
                                    trapezoidColor);
  write("Volumes_TrapezoidVolume");

  return cStream.str();
}

}  // namespace VolumeVisualization
}  // namespace Acts
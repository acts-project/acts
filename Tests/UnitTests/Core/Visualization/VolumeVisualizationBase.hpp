// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Visualization/GeometryVisualization.hpp"
#include "Acts/Visualization/IVisualization.hpp"

#include "Acts/Geometry/AbstractVolume.hpp"
#include "Acts/Geometry/CuboidVolumeBounds.hpp"
#include "Acts/Geometry/CylinderVolumeBounds.hpp"
#include "Acts/Geometry/GenericCuboidVolumeBounds.hpp"
#include "Acts/Geometry/TrapezoidVolumeBounds.hpp"

#include <fstream>
#include <string>

namespace Acts {
namespace VolumeVisualization {
/// Helper method to visualiza all types of surfaces
///
/// @param helper The visualziation helper
/// @param triangulate The directive whether to create triangular meshes
/// @param tag The test tag (mode) identification
/// @param suffix The file suffix for writing
/// @param msuffix the (optional) material file suffix
static inline void test(IVisualization& helper, bool triangulate,
                        const std::string& tag) {
  auto gctx = GeometryContext();
  auto identity = std::make_shared<Transform3D>(Transform3D::Identity());
  std::ofstream stream;

  double halfPhiSector = M_PI / 4.;

  /// Helper method to prepare the streams & helpers
  /// @param path is the file path
  /// @param clear ist he indicator to clear the helper
  auto write = [&](const std::string& path, bool clear = true) -> void {
    std::string wpath = path + tag;
    helper.write(wpath);
    if (clear) {
      helper.clear();
    }
  };

  // ---------------------------------------------------
  // Cuboid surface section
  IVisualization::ColorType boxColor = {0, 0, 255};
  ;
  auto box = std::make_shared<CuboidVolumeBounds>(4., 3., 6.);
  auto cuboid = std::make_shared<AbstractVolume>(identity, box);
  Visualization::drawVolume(helper, *cuboid, gctx, Transform3D::Identity(), 72,
                            triangulate, boxColor);
  write("CuboidVolume");

  //----------------------------------------------------
  // Cylinder volume section
  IVisualization::ColorType cylinderColor = {0, 196, 252};

  double cylinderInnerR = 1.;
  double cylinderOuterR = 5.;
  double cylinderHalfZ = 10.;

  auto fullCylinder =
      std::make_shared<CylinderVolumeBounds>(0., cylinderOuterR, cylinderHalfZ);
  auto cylinder = std::make_shared<AbstractVolume>(identity, fullCylinder);
  Visualization::drawVolume(helper, *cylinder, gctx, Transform3D::Identity(),
                            72, triangulate, cylinderColor);
  write("CylinderVolumeFull");

  auto tubeCylinder = std::make_shared<CylinderVolumeBounds>(
      cylinderInnerR, cylinderOuterR, cylinderHalfZ);
  cylinder = std::make_shared<AbstractVolume>(identity, tubeCylinder);
  Visualization::drawVolume(helper, *cylinder, gctx, Transform3D::Identity(),
                            72, triangulate, cylinderColor);
  write("CylinderVolumeTube");

  tubeCylinder = std::make_shared<CylinderVolumeBounds>(
      cylinderInnerR, cylinderOuterR, cylinderHalfZ, halfPhiSector);
  cylinder = std::make_shared<AbstractVolume>(identity, tubeCylinder);
  Visualization::drawVolume(helper, *cylinder, gctx, Transform3D::Identity(),
                            72, triangulate, cylinderColor);
  write("CylinderVolumeTubeSector");

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
  Visualization::drawVolume(helper, *generic, gctx, Transform3D::Identity(), 72,
                            triangulate, genericColor);
  write("GenericCuboidVolume");

  //----------------------------------------------------
  // Trapezoid volume section
  IVisualization::ColorType trapezoidColor = {23, 196, 22};
  auto trapezoid = std::make_shared<TrapezoidVolumeBounds>(2., 4., 5., 6.);
  auto trapezoidVolume = std::make_shared<AbstractVolume>(identity, trapezoid);
  Visualization::drawVolume(helper, *trapezoidVolume, gctx,
                            Transform3D::Identity(), 72, triangulate,
                            trapezoidColor);
  write("TrapezoidVolume");
}

}  // namespace VolumeVisualization
}  // namespace Acts
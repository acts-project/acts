// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Visualization/GeometryVisualization.hpp"
#include "Acts/Visualization/IVisualization.hpp"

#include "Acts/Surfaces/ConeBounds.hpp"
#include "Acts/Surfaces/ConeSurface.hpp"
#include "Acts/Surfaces/CylinderBounds.hpp"
#include "Acts/Surfaces/CylinderSurface.hpp"
#include "Acts/Surfaces/DiamondBounds.hpp"
#include "Acts/Surfaces/DiscSurface.hpp"
#include "Acts/Surfaces/DiscTrapezoidBounds.hpp"
#include "Acts/Surfaces/EllipseBounds.hpp"
#include "Acts/Surfaces/LineBounds.hpp"
#include "Acts/Surfaces/RadialBounds.hpp"
#include "Acts/Surfaces/RectangleBounds.hpp"
#include "Acts/Surfaces/StrawSurface.hpp"

#include "Acts/Geometry/GeometryContext.hpp"

#include <fstream>
#include <string>

namespace Acts {
namespace SurfaceVisualization {

/// Helper method to visualiza all types of surfaces
///
/// @param helper The visualziation helper
/// @param triangulate The directive whether to create triangular meshes
/// @param suffix The file suffix for writing
static inline void test(IVisualization& helper, bool triangulate,
                        const std::string& suffix) {
  auto gctx = GeometryContext();
  auto identity = std::make_shared<Transform3D>(Transform3D::Identity());
  std::ofstream stream;

  double halfPhiSector = M_PI / 4.;
  double centralPhi = M_PI / 2.;

  //----------------------------------------------------
  // Cone Surface section
  IVisualization::ColorType coneColor = {252, 160, 0};

  double coneAlpha = 0.245;
  double coneMinZ = 0.;
  double coneMaxZ = 10.;
  // Full Cone
  stream.open(std::string("ConeSurface") + suffix);
  auto coneBounds = std::make_shared<ConeBounds>(coneAlpha, coneMinZ, coneMaxZ);
  auto cone = Surface::makeShared<ConeSurface>(identity, coneBounds);
  Visualization::drawSurface(helper, *cone, gctx, Transform3D::Identity(), 72,
                             triangulate, coneColor);
  helper.write(stream);
  stream.close();
  // Sectoral Cone
  helper.clear();
  stream.open(std::string("ConeSurfaceSector") + suffix);
  coneBounds = std::make_shared<ConeBounds>(coneAlpha, coneMinZ, coneMaxZ,
                                            halfPhiSector);
  cone = Surface::makeShared<ConeSurface>(identity, coneBounds);
  Visualization::drawSurface(helper, *cone, gctx, Transform3D::Identity(), 72,
                             triangulate, coneColor);
  helper.write(stream);
  stream.close();
  // Sectoral Cone Shifted
  helper.clear();
  stream.open(std::string("ConeSurfaceSectorShifted") + suffix);
  coneBounds = std::make_shared<ConeBounds>(coneAlpha, coneMinZ, coneMaxZ,
                                            halfPhiSector, centralPhi);
  cone = Surface::makeShared<ConeSurface>(identity, coneBounds);
  Visualization::drawSurface(helper, *cone, gctx, Transform3D::Identity(), 72,
                             triangulate, coneColor);
  helper.write(stream);
  stream.close();

  //----------------------------------------------------
  // Cylinder surface section
  IVisualization::ColorType cylinderColor = {0, 196, 252};

  double cylinderRadius = 5.;
  double cylinderHalfZ = 10.;

  // Full Cylinder
  helper.clear();
  stream.open(std::string("CylinderSurface") + suffix);
  auto cylinderBounds =
      std::make_shared<CylinderBounds>(cylinderRadius, cylinderHalfZ);
  auto cylinder =
      Surface::makeShared<CylinderSurface>(identity, cylinderBounds);
  Visualization::drawSurface(helper, *cylinder, gctx, Transform3D::Identity(),
                             72, triangulate, cylinderColor);
  helper.write(stream);
  stream.close();
  // Sectoral Cone
  helper.clear();
  stream.open(std::string("CylinderSurfaceSector") + suffix);
  cylinderBounds = std::make_shared<CylinderBounds>(
      cylinderRadius, cylinderHalfZ, halfPhiSector);
  cylinder = Surface::makeShared<CylinderSurface>(identity, cylinderBounds);
  Visualization::drawSurface(helper, *cylinder, gctx, Transform3D::Identity(),
                             72, triangulate, cylinderColor);
  helper.write(stream);
  stream.close();
  // Sectoral Cone Shifted
  helper.clear();
  stream.open(std::string("CylinderSurfaceSectorShifted") + suffix);
  cylinderBounds = std::make_shared<CylinderBounds>(
      cylinderRadius, cylinderHalfZ, halfPhiSector, centralPhi);
  cylinder = Surface::makeShared<CylinderSurface>(identity, cylinderBounds);
  Visualization::drawSurface(helper, *cylinder, gctx, Transform3D::Identity(),
                             72, triangulate, cylinderColor);
  helper.write(stream);
  stream.close();

  //----------------------------------------------------
  // Disc Surface section
  IVisualization::ColorType discColor = {126, 252, 0};

  double discRmin = 5.;
  double discRmax = 10.;

  std::vector<std::shared_ptr<DiscSurface> > radialSurfaces;

  // Full Disc
  helper.clear();
  stream.open(std::string("DiscSurfaceFull") + suffix);
  auto radialBounds = std::make_shared<RadialBounds>(0., discRmax);
  auto disc = Surface::makeShared<DiscSurface>(identity, radialBounds);
  radialSurfaces.push_back(disc);
  Visualization::drawSurface(helper, *disc, gctx, Transform3D::Identity(), 72,
                             triangulate, discColor);
  helper.write(stream);
  stream.close();

  // Full Sectoral Disc
  helper.clear();
  stream.open(std::string("DiscSurfaceFullSector") + suffix);
  radialBounds = std::make_shared<RadialBounds>(0., discRmax, halfPhiSector);
  disc = Surface::makeShared<DiscSurface>(identity, radialBounds);
  radialSurfaces.push_back(disc);
  Visualization::drawSurface(helper, *disc, gctx, Transform3D::Identity(), 72,
                             triangulate, discColor);
  helper.write(stream);
  stream.close();

  // Full Sectoral Shifted Disc
  helper.clear();
  stream.open(std::string("DiscSurfaceFullSectorShifted") + suffix);
  radialBounds =
      std::make_shared<RadialBounds>(0., discRmax, halfPhiSector, centralPhi);
  disc = Surface::makeShared<DiscSurface>(identity, radialBounds);
  radialSurfaces.push_back(disc);
  Visualization::drawSurface(helper, *disc, gctx, Transform3D::Identity(), 72,
                             triangulate, discColor);
  helper.write(stream);
  stream.close();

  // Full Ring
  helper.clear();
  stream.open(std::string("DiscSurfaceRing") + suffix);
  radialBounds = std::make_shared<RadialBounds>(discRmin, discRmax);
  disc = Surface::makeShared<DiscSurface>(identity, radialBounds);
  radialSurfaces.push_back(disc);
  Visualization::drawSurface(helper, *disc, gctx, Transform3D::Identity(), 72,
                             triangulate, discColor);
  helper.write(stream);
  stream.close();

  // Full Sectoral Rin g
  helper.clear();
  stream.open(std::string("DiscSurfaceRingSector") + suffix);
  radialBounds =
      std::make_shared<RadialBounds>(discRmin, discRmax, halfPhiSector);
  disc = Surface::makeShared<DiscSurface>(identity, radialBounds);
  radialSurfaces.push_back(disc);
  Visualization::drawSurface(helper, *disc, gctx, Transform3D::Identity(), 72,
                             triangulate, discColor);
  helper.write(stream);
  stream.close();

  // Full Sectoral Shifted Ring
  helper.clear();
  stream.open(std::string("DiscSurfaceRingSectorShifted") + suffix);
  radialBounds = std::make_shared<RadialBounds>(discRmin, discRmax,
                                                halfPhiSector, centralPhi);
  disc = Surface::makeShared<DiscSurface>(identity, radialBounds);
  radialSurfaces.push_back(disc);
  Visualization::drawSurface(helper, *disc, gctx, Transform3D::Identity(), 72,
                             triangulate, discColor);
  helper.write(stream);
  stream.close();

  // All in one
  helper.clear();
  std::vector<Transform3D> sixDiscs = {
      Transform3D(Translation3D{-2.5 * discRmax, 1.5 * discRmax, 0.}),
      Transform3D(Translation3D{0., 1.5 * discRmax, 0.}),
      Transform3D(Translation3D{2.5 * discRmax, 1.5 * discRmax, 0.}),
      Transform3D(Translation3D{-2.5 * discRmax, -1.5 * discRmax, 0.}),
      Transform3D(Translation3D{0., -1.5 * discRmax, 0.}),
      Transform3D(Translation3D{2.5 * discRmax, -1.5 * discRmax, 0.})};
  stream.open(std::string("DiscSurfacesOptions") + suffix);
  for (size_t ir = 0; ir < radialSurfaces.size(); ++ir) {
    Visualization::drawSurface(helper, *radialSurfaces[ir], gctx, sixDiscs[ir],
                               72, triangulate, discColor);
  }
  helper.write(stream);
  stream.close();
}

}  // namespace SurfaceVisualization
}  // namespace Acts

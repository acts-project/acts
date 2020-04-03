// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Visualization/GeometryVisualization.hpp"
#include "Acts/Visualization/IVisualization.hpp"

#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Surfaces/AnnulusBounds.hpp"
#include "Acts/Surfaces/ConeBounds.hpp"
#include "Acts/Surfaces/ConeSurface.hpp"
#include "Acts/Surfaces/ConvexPolygonBounds.hpp"
#include "Acts/Surfaces/CylinderBounds.hpp"
#include "Acts/Surfaces/CylinderSurface.hpp"
#include "Acts/Surfaces/DiamondBounds.hpp"
#include "Acts/Surfaces/DiscSurface.hpp"
#include "Acts/Surfaces/DiscTrapezoidBounds.hpp"
#include "Acts/Surfaces/EllipseBounds.hpp"
#include "Acts/Surfaces/LineBounds.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/Surfaces/RadialBounds.hpp"
#include "Acts/Surfaces/RectangleBounds.hpp"
#include "Acts/Surfaces/StrawSurface.hpp"
#include "Acts/Surfaces/TrapezoidBounds.hpp"

#include <fstream>
#include <string>

namespace Acts {
namespace SurfaceVisualization {

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
  double centralPhi = M_PI / 2.;

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

  //----------------------------------------------------
  // Cone Surface section
  IVisualization::ColorType coneColor = {252, 160, 0};
  std::vector<std::shared_ptr<ConeSurface>> coneSurfaces;

  double coneAlpha = 0.245;
  double coneMinZ = 0.;
  double coneCutZ = 2.;
  double coneMaxZ = 10.;
  // Full Cone
  auto coneBounds =
      std::make_shared<ConeBounds>(coneAlpha, -coneCutZ, coneMaxZ);
  auto cone = Surface::makeShared<ConeSurface>(identity, coneBounds);
  coneSurfaces.push_back(cone);
  Visualization::drawSurface(helper, *cone, gctx, Transform3D::Identity(), 72,
                             triangulate, coneColor);
  write("ConeSurface");

  // Sectoral Cone
  coneBounds = std::make_shared<ConeBounds>(coneAlpha, coneMinZ, coneMaxZ,
                                            halfPhiSector);
  cone = Surface::makeShared<ConeSurface>(identity, coneBounds);
  coneSurfaces.push_back(cone);
  Visualization::drawSurface(helper, *cone, gctx, Transform3D::Identity(), 72,
                             triangulate, coneColor);
  write("ConeSurfaceSector");

  // Sectoral Cone Shifted
  coneBounds = std::make_shared<ConeBounds>(coneAlpha, coneCutZ, coneMaxZ,
                                            halfPhiSector, centralPhi);
  cone = Surface::makeShared<ConeSurface>(identity, coneBounds);
  coneSurfaces.push_back(cone);
  Visualization::drawSurface(helper, *cone, gctx, Transform3D::Identity(), 72,
                             triangulate, coneColor);
  write("ConeSurfaceSectorShifted");

  // All in one for radial bounds
  std::vector<Transform3D> threeCones = {
      Transform3D(Translation3D{-0.5 * coneMaxZ, 0., 0.}),
      Transform3D(Translation3D{0., 0., 0.}),
      Transform3D(Translation3D{0.75 * coneMaxZ, 0., 0.})};

  for (size_t ic = 0; ic < coneSurfaces.size(); ++ic) {
    Visualization::drawSurface(helper, *coneSurfaces[ic], gctx, threeCones[ic],
                               72, triangulate, coneColor);
  }
  write("All_ConeSurfaces");

  //----------------------------------------------------
  // Cylinder surface section
  IVisualization::ColorType cylinderColor = {0, 196, 252};
  std::vector<std::shared_ptr<CylinderSurface>> cylinderSurfaces;

  double cylinderRadius = 5.;
  double cylinderHalfZ = 10.;

  // Full Cylinder
  auto cylinderBounds =
      std::make_shared<CylinderBounds>(cylinderRadius, cylinderHalfZ);
  auto cylinder =
      Surface::makeShared<CylinderSurface>(identity, cylinderBounds);
  cylinderSurfaces.push_back(cylinder);
  Visualization::drawSurface(helper, *cylinder, gctx, Transform3D::Identity(),
                             72, triangulate, cylinderColor);
  write("CylinderSurface");
  // Sectoral Cone
  cylinderBounds = std::make_shared<CylinderBounds>(
      cylinderRadius, cylinderHalfZ, halfPhiSector);
  cylinder = Surface::makeShared<CylinderSurface>(identity, cylinderBounds);
  cylinderSurfaces.push_back(cylinder);
  Visualization::drawSurface(helper, *cylinder, gctx, Transform3D::Identity(),
                             72, triangulate, cylinderColor);
  write("CylinderSurfaceSector");

  // Sectoral Cone Shifted
  cylinderBounds = std::make_shared<CylinderBounds>(
      cylinderRadius, cylinderHalfZ, halfPhiSector, centralPhi);
  cylinder = Surface::makeShared<CylinderSurface>(identity, cylinderBounds);
  cylinderSurfaces.push_back(cylinder);
  Visualization::drawSurface(helper, *cylinder, gctx, Transform3D::Identity(),
                             72, triangulate, cylinderColor);
  write("CylinderSurfaceSectorShifted");

  // All in one for radial bounds
  std::vector<Transform3D> threeCylinders = {
      Transform3D(Translation3D{-2 * cylinderRadius, 0., 0.}),
      Transform3D(Translation3D{0., 0., 0.}),
      Transform3D(Translation3D{2.5 * cylinderRadius, 0., 0.})};

  for (size_t ic = 0; ic < cylinderSurfaces.size(); ++ic) {
    Visualization::drawSurface(helper, *cylinderSurfaces[ic], gctx,
                               threeCylinders[ic], 72, triangulate,
                               cylinderColor);
  }
  write("All_CylinderSurfaces");

  /// ------------- planar bounding box
  /// @param name of the file
  auto writeBoundingBox2D = [&](const RectangleBounds& rBounds,
                                const std::string& path) -> void {
    IVisualization::ColorType bbColor = {126, 126, 126};

    std::string bbPath = path + tag + "_bbox";

    auto bbBounds = std::make_shared<RectangleBounds>(rBounds);
    auto bbSurface = Surface::makeShared<PlaneSurface>(identity, bbBounds);
    Visualization::drawSurface(helper, *bbSurface, gctx,
                               Transform3D::Identity(), 72, triangulate,
                               bbColor);
    write(bbPath);
  };

  //----------------------------------------------------
  // Disc Surface section
  IVisualization::ColorType discColor = {126, 252, 0};

  double discRmin = 5.;
  double discRmax = 10.;

  std::vector<std::shared_ptr<DiscSurface>> radialSurfaces;

  // Full Disc
  auto radialBounds = std::make_shared<RadialBounds>(0., discRmax);
  auto disc = Surface::makeShared<DiscSurface>(identity, radialBounds);
  radialSurfaces.push_back(disc);
  Visualization::drawSurface(helper, *disc, gctx, Transform3D::Identity(), 72,
                             triangulate, discColor);
  write("DiscSurfaceFull");

  // Full Sectoral Disc
  radialBounds = std::make_shared<RadialBounds>(0., discRmax, halfPhiSector);
  disc = Surface::makeShared<DiscSurface>(identity, radialBounds);
  radialSurfaces.push_back(disc);
  Visualization::drawSurface(helper, *disc, gctx, Transform3D::Identity(), 72,
                             triangulate, discColor);
  write("DiscSurfaceFullSector");

  // Full Sectoral Shifted Disc
  radialBounds =
      std::make_shared<RadialBounds>(0., discRmax, halfPhiSector, centralPhi);
  disc = Surface::makeShared<DiscSurface>(identity, radialBounds);
  radialSurfaces.push_back(disc);
  Visualization::drawSurface(helper, *disc, gctx, Transform3D::Identity(), 72,
                             triangulate, discColor);
  write("DiscSurfaceFullSectorShifted");

  // Full Ring
  radialBounds = std::make_shared<RadialBounds>(discRmin, discRmax);
  disc = Surface::makeShared<DiscSurface>(identity, radialBounds);
  radialSurfaces.push_back(disc);
  Visualization::drawSurface(helper, *disc, gctx, Transform3D::Identity(), 72,
                             triangulate, discColor);
  write("DiscSurfaceRing");

  // Full Sectoral Rin g
  radialBounds =
      std::make_shared<RadialBounds>(discRmin, discRmax, halfPhiSector);
  disc = Surface::makeShared<DiscSurface>(identity, radialBounds);
  radialSurfaces.push_back(disc);
  Visualization::drawSurface(helper, *disc, gctx, Transform3D::Identity(), 72,
                             triangulate, discColor);
  write("DiscSurfaceRingSector");

  // Full Sectoral Shifted Ring
  radialBounds = std::make_shared<RadialBounds>(discRmin, discRmax,
                                                halfPhiSector, centralPhi);
  disc = Surface::makeShared<DiscSurface>(identity, radialBounds);
  radialSurfaces.push_back(disc);
  Visualization::drawSurface(helper, *disc, gctx, Transform3D::Identity(), 72,
                             triangulate, discColor);
  write("DiscSurfaceRingSectorShifted");

  // All in one for radial bounds
  std::vector<Transform3D> sixDiscs = {
      Transform3D(Translation3D{-2.0 * discRmax, 1.5 * discRmax, 0.}),
      Transform3D(Translation3D{0., 1.5 * discRmax, 0.}),
      Transform3D(Translation3D{2.5 * discRmax, 1.5 * discRmax, 0.}),
      Transform3D(Translation3D{-2.0 * discRmax, -1.5 * discRmax, 0.}),
      Transform3D(Translation3D{0., -1.5 * discRmax, 0.}),
      Transform3D(Translation3D{2.5 * discRmax, -1.5 * discRmax, 0.})};
  for (size_t ir = 0; ir < radialSurfaces.size(); ++ir) {
    Visualization::drawSurface(helper, *radialSurfaces[ir], gctx, sixDiscs[ir],
                               72, triangulate, discColor);
  }
  write("All_DiscSurfaces_RadialBounds");

  std::vector<std::shared_ptr<DiscSurface>> anomalDiscSurfaces;
  IVisualization::ColorType discAnomalColor = {94, 186, 2};

  double annulusMinPhi = 0.75;
  double annulusMaxPhi = 1.35;
  Vector2D offset(-4., 2.);
  auto annulus = std::make_shared<AnnulusBounds>(
      discRmin, discRmax, annulusMinPhi, annulusMaxPhi, offset);
  disc = Surface::makeShared<DiscSurface>(identity, annulus);
  anomalDiscSurfaces.push_back(disc);
  Visualization::drawSurface(helper, *disc, gctx, Transform3D::Identity(), 72,
                             triangulate, discAnomalColor);
  write("DiscAnulusBounds");

  double discTrapezoidHxRmin = 3.;
  double discTrapezoidHxRmax = 6.;
  auto discTrapezoid = std::make_shared<DiscTrapezoidBounds>(
      discTrapezoidHxRmin, discTrapezoidHxRmax, discRmin, discRmax);
  disc = Surface::makeShared<DiscSurface>(identity, discTrapezoid);
  anomalDiscSurfaces.push_back(disc);
  Visualization::drawSurface(helper, *disc, gctx, Transform3D::Identity(), 72,
                             triangulate, discAnomalColor);
  write("DiscTrapezoidBounds");

  // All in one for radial bounds
  std::vector<Transform3D> twoAnomalDiscs = {
      Transform3D(Translation3D{-5., 0., 0.}),
      Transform3D(Translation3D{5., 0., 0.})};
  for (size_t id = 0; id < anomalDiscSurfaces.size(); ++id) {
    Visualization::drawSurface(helper, *anomalDiscSurfaces[id], gctx,
                               sixDiscs[id], 72, triangulate, discAnomalColor);
  }
  write("All_DiscSurfaces_AnomalBounds");

  //----------------------------------------------------
  // Plane Surface section
  IVisualization::ColorType planeColor = {0, 0, 255};
  std::vector<std::shared_ptr<PlaneSurface>> planarSurfaces;

  // Ellipse shaped : Full Ellipse
  double ellipseR0min = 2;
  double ellipseR0max = 4;
  double ellipseR1min = 3;
  double ellipseR1max = 6;
  std::string name = "PlaneSurfaceEllipse";
  auto ellipse =
      std::make_shared<EllipseBounds>(0., 0., ellipseR1min, ellipseR1max);
  auto plane = Surface::makeShared<PlaneSurface>(identity, ellipse);
  planarSurfaces.push_back(plane);
  Visualization::drawSurface(helper, *plane, gctx, Transform3D::Identity(), 72,
                             triangulate, planeColor);
  write(name);
  writeBoundingBox2D(ellipse->boundingBox(), name);

  // Ellipse shaped : Ring Ellipse
  name = "PlaneSurfaceEllipseRing";
  ellipse = std::make_shared<EllipseBounds>(ellipseR0min, ellipseR0max,
                                            ellipseR1min, ellipseR1max);
  plane = Surface::makeShared<PlaneSurface>(identity, ellipse);
  planarSurfaces.push_back(plane);
  Visualization::drawSurface(helper, *plane, gctx, Transform3D::Identity(), 72,
                             triangulate, planeColor);
  write(name);
  writeBoundingBox2D(ellipse->boundingBox(), name);

  // Ellipse shaped : Ring Ellipse Sector
  name = "PlaneSurfaceEllipseRingSector";
  ellipse = std::make_shared<EllipseBounds>(
      ellipseR0min, ellipseR0max, ellipseR1min, ellipseR1max, halfPhiSector);
  plane = Surface::makeShared<PlaneSurface>(identity, ellipse);
  planarSurfaces.push_back(plane);
  Visualization::drawSurface(helper, *plane, gctx, Transform3D::Identity(), 72,
                             triangulate, planeColor);
  write(name);
  writeBoundingBox2D(ellipse->boundingBox(), name);

  // ConvexPolygon shaped example: Triangle
  name = "PlaneSurfaceTriangleRegular";
  std::vector<Vector2D> tvertices = {{-3, -1.5}, {3, -1.5}, {0, 4.5}};
  auto triangle = std::make_shared<ConvexPolygonBounds<3>>(tvertices);
  plane = Surface::makeShared<PlaneSurface>(identity, triangle);
  planarSurfaces.push_back(plane);
  Visualization::drawSurface(helper, *plane, gctx, Transform3D::Identity(), 72,
                             triangulate, planeColor);
  write(name);
  writeBoundingBox2D(triangle->boundingBox(), name);

  // ConvexPolygon shaped example: Triangle
  name = "PlaneSurfaceTriangleGeneral";
  tvertices = {{-1., 4.5}, {4, 6.5}, {3, 8.5}};
  triangle = std::make_shared<ConvexPolygonBounds<3>>(tvertices);
  plane = Surface::makeShared<PlaneSurface>(identity, triangle);
  planarSurfaces.push_back(plane);
  Visualization::drawSurface(helper, *plane, gctx, Transform3D::Identity(), 72,
                             triangulate, planeColor);
  write(name);
  writeBoundingBox2D(triangle->boundingBox(), name);

  // ConvexPolygon shaped example: Triangle
  name = "PlaneSurfaceConvexPolygonGeneral";
  tvertices = {{-1., 4.5}, {4, 6.5}, {6, 8.5}, {0, 10.5}, {-3, 6.2}};
  auto dynamicpolygon =
      std::make_shared<ConvexPolygonBounds<PolygonDynamic>>(tvertices);
  plane = Surface::makeShared<PlaneSurface>(identity, dynamicpolygon);
  planarSurfaces.push_back(plane);
  Visualization::drawSurface(helper, *plane, gctx, Transform3D::Identity(), 72,
                             triangulate, planeColor);
  write(name);
  writeBoundingBox2D(dynamicpolygon->boundingBox(), name);

  // Diamond shaped
  name = "PlaneSurfaceDiamond";
  auto diamond = std::make_shared<DiamondBounds>(3., 6., 2., 2., 4.);
  plane = Surface::makeShared<PlaneSurface>(identity, diamond);
  planarSurfaces.push_back(plane);
  Visualization::drawSurface(helper, *plane, gctx, Transform3D::Identity(), 72,
                             triangulate, planeColor);
  write(name);
  writeBoundingBox2D(diamond->boundingBox(), name);

  // Rectangle plane
  name = "PlaneSurfaceRectangle";
  auto rectangle = std::make_shared<RectangleBounds>(2., 3.);
  plane = Surface::makeShared<PlaneSurface>(identity, rectangle);
  planarSurfaces.push_back(plane);
  Visualization::drawSurface(helper, *plane, gctx, Transform3D::Identity(), 72,
                             triangulate, planeColor);
  write(name);
  writeBoundingBox2D(rectangle->boundingBox(), name);

  // Off-centered Rectangle plane:
  name = "PlaneSurfaceRectangleOffcentered";
  rectangle =
      std::make_shared<RectangleBounds>(Vector2D{1., 2.}, Vector2D{15., 12.});
  plane = Surface::makeShared<PlaneSurface>(identity, rectangle);
  Visualization::drawSurface(helper, *plane, gctx, Transform3D::Identity(), 72,
                             triangulate, planeColor);
  write(name);
  writeBoundingBox2D(rectangle->boundingBox(), name);

  // Off-centered Rectangle plane:
  name = "PlaneSurfaceTrapezoid";
  auto trapezoid = std::make_shared<TrapezoidBounds>(2., 5., 3.);
  plane = Surface::makeShared<PlaneSurface>(identity, trapezoid);
  planarSurfaces.push_back(plane);
  Visualization::drawSurface(helper, *plane, gctx, Transform3D::Identity(), 72,
                             triangulate, planeColor);
  write(name);
  writeBoundingBox2D(trapezoid->boundingBox(), name);

  // All planes
  std::vector<Transform3D> ninePlanes = {
      Transform3D(Translation3D{-10., -10., 0.}),
      Transform3D(Translation3D{0., -10., 0.}),
      Transform3D(Translation3D{10., -10., 0.}),
      Transform3D(Translation3D{-10., 0., 0.}),
      Transform3D(Translation3D{0., -6., 0.}),
      Transform3D(Translation3D{10., -8., 0.}),
      Transform3D(Translation3D{-10, 10., 0.}),
      Transform3D(Translation3D{0., 10., 0.}),
      Transform3D(Translation3D{10., 10., 0.})};
  for (size_t ip = 0; ip < planarSurfaces.size(); ++ip) {
    Visualization::drawSurface(helper, *planarSurfaces[ip], gctx,
                               ninePlanes[ip], 72, triangulate, planeColor);
  }
  write("All_PlaneSurfaces");

  //----------------------------------------------------
  // Straw Surface section
  IVisualization::ColorType strawColor = {255, 0, 0};

  name = "StrawSurface";
  auto tube = std::make_shared<LineBounds>(2., 20.);
  auto straw = Surface::makeShared<StrawSurface>(identity, tube);
  Visualization::drawSurface(helper, *straw, gctx, Transform3D::Identity(), 72,
                             triangulate, strawColor);
  write(name);
}

}  // namespace SurfaceVisualization
}  // namespace Acts

// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Surfaces/AnnulusBounds.hpp"
#include "Acts/Surfaces/ConvexPolygonBounds.hpp"
#include "Acts/Surfaces/DiamondBounds.hpp"
#include "Acts/Surfaces/DiscSurface.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/Surfaces/RadialBounds.hpp"
#include "Acts/Surfaces/RectangleBounds.hpp"
#include "Acts/Surfaces/TrapezoidBounds.hpp"
#include "ActsPlugins/ActSVG/SurfaceSvgConverter.hpp"
#include "ActsPlugins/ActSVG/SvgUtils.hpp"

#include <fstream>
#include <numbers>

using namespace Acts;
using namespace ActsPlugins;

namespace ActsTests {

/// Helper to run planar Test
///
/// @param surface
/// @param identification
void runPlanarTests(const Surface& surface, const Svg::Style& style,
                    const std::string& identification) {
  // Default Geometry context
  auto geoCtx = GeometryContext::dangerouslyDefaultConstruct();

  using SurfaceOptions = Svg::SurfaceConverter::Options;

  SurfaceOptions sOptions;
  sOptions.style = style;
  sOptions.templateSurface = true;

  // Svg proto object & actual object
  auto svgTemplate = Svg::SurfaceConverter::convert(geoCtx, surface, sOptions);
  auto xyTemplate = Svg::View::xy(svgTemplate, identification + "_template");
  Svg::toFile({xyTemplate}, xyTemplate._id + ".svg");
  // Positioned
  sOptions.templateSurface = false;
  auto svgObject = Svg::SurfaceConverter::convert(geoCtx, surface, sOptions);
  auto xyObject = Svg::View::xy(svgObject, identification);
  auto xyAxes = Svg::axesXY(static_cast<double>(xyObject._x_range[0]),
                            static_cast<double>(xyObject._x_range[1]),
                            static_cast<double>(xyObject._y_range[0]),
                            static_cast<double>(xyObject._y_range[1]));

  Svg::toFile({xyObject, xyAxes}, xyObject._id + ".svg");
  // As sheet
  auto svgSheet = Svg::Sheet::xy(svgTemplate, identification + "_sheet");
  Svg::toFile({svgSheet}, svgSheet._id + ".svg");
}

BOOST_AUTO_TEST_SUITE(ActSvg)

BOOST_AUTO_TEST_CASE(PlanarSurfaces) {
  // Planar style
  Svg::Style planarStyle;
  planarStyle.fillColor = {51, 153, 255};
  planarStyle.fillOpacity = 0.75;
  planarStyle.highlightColor = {255, 153, 51};
  planarStyle.highlights = {"mouseover", "mouseout"};
  planarStyle.strokeWidth = 0.5;
  planarStyle.quarterSegments = 0u;

  // Rectangle case
  auto rectangleBounds = std::make_shared<RectangleBounds>(36., 64.);
  auto transform = Transform3::Identity();
  transform.pretranslate(Vector3{20., 20., 100.});
  auto rectanglePlane =
      Surface::makeShared<PlaneSurface>(transform, rectangleBounds);
  runPlanarTests(*rectanglePlane, planarStyle, "rectangle");

  // Trapezoid case:
  auto trapezoidBounds = std::make_shared<TrapezoidBounds>(36., 64., 105.);
  auto trapeozidPlane =
      Surface::makeShared<PlaneSurface>(transform, trapezoidBounds);
  runPlanarTests(*trapeozidPlane, planarStyle, "trapezoid");

  // Trapezoid case shifted and rotated
  double phi = std::numbers::pi / 8.;
  double radius = 150.;
  Vector3 center(radius * std::cos(phi), radius * std::sin(phi), 0.);

  Vector3 localY(std::cos(phi), std::sin(phi), 0.);
  Vector3 localZ(0., 0., 1.);
  Vector3 localX = localY.cross(localZ);
  RotationMatrix3 rotation;
  rotation.col(0) = localX;
  rotation.col(1) = localY;
  rotation.col(2) = localZ;
  transform = Transform3(Translation3(center) * rotation);
  // Create the module surface
  auto trapeozidPlaneTransformed =
      Surface::makeShared<PlaneSurface>(transform, trapezoidBounds);

  runPlanarTests(*trapeozidPlaneTransformed, planarStyle, "trapezoid_rotated");
  // A reference test for the rotated one
  auto geoCtx = GeometryContext::dangerouslyDefaultConstruct();
  actsvg::proto::surface<std::vector<Vector3>> reference;
  reference._vertices =
      trapeozidPlaneTransformed->polyhedronRepresentation(geoCtx, 1u).vertices;
  auto referenceTrapezoid = Svg::View::xy(reference, "trapezoid");
  auto referenceAxes = Svg::axesXY(-200., 200., -200., 200);
  Svg::toFile({referenceTrapezoid, referenceAxes}, "trapezoid_reference.svg");

  // Let's create one with a flipped z-axis
  Vector3 flocalZ(0., 0., -1.);
  Vector3 flocalX = localY.cross(flocalZ);
  RotationMatrix3 frotation;
  frotation.col(0) = flocalX;
  frotation.col(1) = localY;
  frotation.col(2) = flocalZ;
  auto ftransform = Transform3(Translation3(center) * frotation);
  // Create the module surface
  auto ftrapeozidPlaneTransformed =
      Surface::makeShared<PlaneSurface>(ftransform, trapezoidBounds);

  runPlanarTests(*ftrapeozidPlaneTransformed, planarStyle,
                 "flipped_trapezoid_rotated");
  actsvg::proto::surface<std::vector<Vector3>> freference;
  freference._vertices =
      ftrapeozidPlaneTransformed->polyhedronRepresentation(geoCtx, 1u).vertices;

  auto freferenceTrapezoid = Svg::View::xy(freference, "flipped_trapezoid");
  Svg::toFile({freferenceTrapezoid, referenceAxes},
              "flipped_trapezoid_reference.svg");

  // Diamond
  auto diamondBounds = std::make_shared<DiamondBounds>(36., 64., 14., 40., 30.);
  transform = Transform3::Identity();
  auto diamond = Surface::makeShared<PlaneSurface>(transform, diamondBounds);
  runPlanarTests(*diamond, planarStyle, "diamond");

  // ConvexPolygon
  std::vector<Vector2> vertices = {
      {-10., -10.}, {10., -15.}, {20., 5.}, {-5., 15.}, {-12, 0.}};
  auto polygonBounds = std::make_shared<ConvexPolygonBounds<5u>>(vertices);
  auto polygon = Surface::makeShared<PlaneSurface>(transform, polygonBounds);
  runPlanarTests(*polygon, planarStyle, "polygon");
}

BOOST_AUTO_TEST_CASE(DiscSurfaces) {
  // Planar style
  Svg::Style discStyle;
  discStyle.fillColor = {0, 204, 153};
  discStyle.fillOpacity = 0.75;
  discStyle.highlightColor = {153, 204, 0};
  discStyle.highlights = {"mouseover", "mouseout"};
  discStyle.strokeWidth = 0.5;
  discStyle.quarterSegments = 72u;

  auto transform = Transform3::Identity();
  transform.pretranslate(Vector3{20., 20., 100.});

  // Full disc case
  auto fullDiscBounds = std::make_shared<RadialBounds>(0., 64.);
  auto fullDisc = Surface::makeShared<DiscSurface>(transform, fullDiscBounds);
  runPlanarTests(*fullDisc, discStyle, "full_disc");

  // Full ring case:
  auto fullRingBounds = std::make_shared<RadialBounds>(36., 64.);
  auto fullRing = Surface::makeShared<DiscSurface>(transform, fullRingBounds);
  runPlanarTests(*fullRing, discStyle, "full_ring");

  // Sectorial disc case
  auto sectoralDiscBounds = std::make_shared<RadialBounds>(
      0., 64., std::numbers::pi / 4., std::numbers::pi / 2.);
  auto sectoralDisc =
      Surface::makeShared<DiscSurface>(transform, sectoralDiscBounds);
  runPlanarTests(*sectoralDisc, discStyle, "full_disc");

  // Annulus shape
  double minRadius = 7.2;
  double maxRadius = 12.0;
  double minPhi = 0.74195;
  double maxPhi = 1.33970;

  Vector2 offset{-3., 2.};

  auto annulusDiscBounds = std::make_shared<AnnulusBounds>(
      minRadius, maxRadius, minPhi, maxPhi, offset);
  auto annulusDisc =
      Surface::makeShared<DiscSurface>(transform, annulusDiscBounds);
  runPlanarTests(*annulusDisc, discStyle, "annulus_disc");
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace ActsTests

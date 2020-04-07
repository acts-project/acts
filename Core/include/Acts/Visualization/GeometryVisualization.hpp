// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Geometry/Polyhedron.hpp"
#include "Acts/Surfaces/ConeBounds.hpp"
#include "Acts/Surfaces/ConeSurface.hpp"
#include "Acts/Surfaces/CylinderBounds.hpp"
#include "Acts/Surfaces/CylinderSurface.hpp"
#include "Acts/Surfaces/DiscSurface.hpp"
#include "Acts/Surfaces/RadialBounds.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Utilities/UnitVectors.hpp"
#include "Acts/Visualization/IVisualization.hpp"

namespace Acts {

namespace Visualization {

/// Helper method to Surface objects
///
/// @param helper [in, out] The visualization helper
/// @param surface The surface to be drawn
/// @param gctx The geometry context for which it is drawn
/// @param transform An option additional transform
/// @param lseg The number of segments for a full arch (if needed)
/// @param triangulate The (optional) boolean diretive to triangulate the faces
/// @param color the (optional) color of the object to be written
static inline void drawSurface(
    IVisualization& helper, const Surface& surface, const GeometryContext& gctx,
    const Transform3D& transform = Transform3D::Identity(), size_t lseg = 72,
    bool triangulate = false,
    const IVisualization::ColorType& color = {120, 120, 120}) {
  // Drawing the polyhedron representation of surfaces
  Polyhedron phedron = surface.polyhedronRepresentation(gctx, lseg);
  phedron.move(transform);
  phedron.draw(helper, triangulate, color);
}

/// Helper method for volume objects
///
/// @tparam volume_t the templated volume class
///
/// @param helper [in, out] The visualization helper
/// @param volume The surface to be drawn
/// @param gctx The geometry context for which it is drawn
/// @param transform An option additional transform
/// @param lseg The number of segments for a full arch (if needed)
/// @param triangulate The (optional) boolean diretive to triangulate the faces
/// @param color the (optional) color of the object to be written
template <typename volume_t>
inline void drawVolume(IVisualization& helper, const volume_t& volume,
                       const GeometryContext& gctx,
                       const Transform3D& transform = Transform3D::Identity(),
                       size_t lseg = 72, bool triangulate = false,
                       const IVisualization::ColorType& color = {120, 120,
                                                                 120}) {
  // Drawing the polyhedron representation of surfaces
  auto bSurfaces = volume.boundarySurfaces();
  for (const auto& bs : bSurfaces) {
    drawSurface(helper, bs->surfaceRepresentation(), gctx, transform, lseg,
                triangulate, color);
  }
}

/// Helper method to draw lines - base for all lines
///
/// @param helper [in, out] The visualization helper
/// @param start The start point
/// @param end The end point
/// @param thickness of the line, if bigger 0, approximated by cylinder
/// @param arrows [ -1 | 0 | 1 | 2 ] = [ start | none | end | both ]
/// @param arrowLength wrt halflength
/// @param arrowWidth wrt thickness
/// @param lseg The number of segments for a full arch (if needed)
/// @param color the (optional) color of the object to be written
static inline void drawSegmentBase(
    IVisualization& helper, const Vector3D& start, const Vector3D& end,
    double thickness, int arrows = 0, double arrowLength = 0.,
    double arrowWidth = 0., size_t lseg = 72,
    const IVisualization::ColorType& color = {120, 120, 120}) {
  // Draw the parameter shaft and cone
  auto direction = Vector3D(end - start).normalized();
  double hlength = 0.5 * Vector3D(end - start).norm();

  auto unitVectors = makeCurvilinearUnitVectors(direction);
  RotationMatrix3D lrotation;
  lrotation.col(0) = unitVectors.first;
  lrotation.col(1) = unitVectors.second;
  lrotation.col(2) = direction;

  Vector3D lcenter = 0.5 * (start + end);
  double alength = arrowLength * (2 * hlength);

  if (arrows == 2) {
    hlength -= alength;
  } else if (arrows != 0) {
    hlength -= 0.5 * alength;
    lcenter -= Vector3D(arrows * 0.5 * alength * direction);
  }

  // Line - draw a line
  if (thickness > 0.) {
    auto ltransform = std::make_shared<Transform3D>(Transform3D::Identity());
    ltransform->prerotate(lrotation);
    ltransform->pretranslate(lcenter);

    auto lbounds = std::make_shared<CylinderBounds>(thickness, hlength);
    auto line = Surface::makeShared<CylinderSurface>(ltransform, lbounds);

    drawSurface(helper, *line, GeometryContext(), Transform3D::Identity(), lseg,
                false, color);
  } else {
    helper.line(start, end, color);
  }

  // Arrowheads - if configured
  if (arrows != 0) {
    double awith = thickness * arrowWidth;
    double alpha = atan2(thickness * arrowWidth, alength);
    auto plateBounds = std::make_shared<RadialBounds>(thickness, awith);

    if (arrows > 0) {
      auto aetransform = std::make_shared<Transform3D>(Transform3D::Identity());
      aetransform->prerotate(lrotation);
      aetransform->pretranslate(end);
      // Arrow cone
      auto coneBounds = std::make_shared<ConeBounds>(alpha, -alength, 0.);
      auto cone = Surface::makeShared<ConeSurface>(aetransform, coneBounds);
      drawSurface(helper, *cone, GeometryContext(), Transform3D::Identity(),
                  lseg, false, color);
      // Arrow end plate
      auto aptransform = std::make_shared<Transform3D>(Transform3D::Identity());
      aptransform->prerotate(lrotation);
      aptransform->pretranslate(Vector3D(end - alength * direction));

      auto plate = Surface::makeShared<DiscSurface>(aptransform, plateBounds);
      drawSurface(helper, *plate, GeometryContext(), Transform3D::Identity(),
                  lseg, false, color);
    }
    if (arrows < 0 or arrows == 2) {
      auto astransform = std::make_shared<Transform3D>(Transform3D::Identity());
      astransform->prerotate(lrotation);
      astransform->pretranslate(start);

      // Arrow cone
      auto coneBounds = std::make_shared<ConeBounds>(alpha, 0., alength);
      auto cone = Surface::makeShared<ConeSurface>(astransform, coneBounds);
      drawSurface(helper, *cone, GeometryContext(), Transform3D::Identity(),
                  lseg, false, color);
      // Arrow end plate
      auto aptransform = std::make_shared<Transform3D>(Transform3D::Identity());
      aptransform->prerotate(lrotation);
      aptransform->pretranslate(Vector3D(start + alength * direction));

      auto plate = Surface::makeShared<DiscSurface>(aptransform, plateBounds);
      drawSurface(helper, *plate, GeometryContext(), Transform3D::Identity(),
                  lseg, false, color);
    }
  }
}

/// Convenience function : line
///
/// @param helper [in, out] The visualization helper
/// @param start The start point
/// @param end The end point
/// @param thickness of the line, if bigger 0, approximated by cylinder
/// @param lseg The number of segments for a full arch (if needed)
/// @param color the (optional) color of the object to be written
static inline void drawSegment(IVisualization& helper, const Vector3D& start,
                               const Vector3D& end, double thickness,
                               size_t lseg = 72,
                               const IVisualization::ColorType& color = {
                                   20, 120, 120}) {
  drawSegmentBase(helper, start, end, thickness, 0, 0., 0., lseg, color);
}

/// Convenience function : arrow pointing back
///
/// @param helper [in, out] The visualization helper
/// @param start The start point
/// @param end The end point
/// @param thickness of the line, if bigger 0, approximated by cylinder
/// @param arrowLength wrt halflength
/// @param arrorWidth wrt thickness
/// @param lseg The number of segments for a full arch (if needed)
/// @param color the (optional) color of the object to be written
static inline void drawArrowBackward(
    IVisualization& helper, const Vector3D& start, const Vector3D& end,
    double thickness, double arrowLength, double arrowWidth, size_t lseg = 72,
    const IVisualization::ColorType& color = {120, 120, 120}) {
  drawSegmentBase(helper, start, end, thickness, -1, arrowLength, arrowWidth,
                  lseg, color);
}

/// Convenience function : arrow pointing forwad
///
/// @param helper [in, out] The visualization helper
/// @param start The start point
/// @param end The end point
/// @param thickness of the line, if bigger 0, approximated by cylinder
/// @param arrowLength wrt halflength
/// @param arrorWidth wrt thickness
/// @param lseg The number of segments for a full arch (if needed)
/// @param color the (optional) color of the object to be written
static inline void drawArrowForward(
    IVisualization& helper, const Vector3D& start, const Vector3D& end,
    double thickness, double arrowLength, double arrowWidth, size_t lseg = 72,
    const IVisualization::ColorType& color = {120, 120, 120}) {
  drawSegmentBase(helper, start, end, thickness, 1, arrowLength, arrowWidth,
                  lseg, color);
}

/// Convenience function : arrow pointing both directions
///
/// @param helper [in, out] The visualization helper
/// @param start The start point
/// @param end The end point
/// @param thickness of the line, if bigger 0, approximated by cylinder
/// @param arrowLength wrt halflength
/// @param arrorWidth wrt thickness
/// @param lseg The number of segments for a full arch (if needed)
/// @param color the (optional) color of the object to be written
static inline void drawArrowsBoth(
    IVisualization& helper, const Vector3D& start, const Vector3D& end,
    double thickness, double arrowLength, double arrowWidth, size_t lseg = 72,
    const IVisualization::ColorType& color = {120, 120, 120}) {
  drawSegmentBase(helper, start, end, thickness, 2, arrowLength, arrowWidth,
                  lseg, color);
}

}  // namespace Visualization

}  // namespace Acts
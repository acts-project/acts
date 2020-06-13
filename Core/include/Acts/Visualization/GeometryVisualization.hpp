// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Geometry/AbstractVolume.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Surfaces/SurfaceArray.hpp"
#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Visualization/IVisualization.hpp"

namespace Acts {

struct GeometryVisualization {
  /// Helper method to Surface objects
  ///
  /// @param helper [in, out] The visualization helper
  /// @param surface The surface to be drawn
  /// @param gctx The geometry context for which it is drawn
  /// @param transform An option additional transform
  /// @param lseg The number of segments for a full arch (if needed)
  /// @param triangulate The (optional) boolean diretive to triangulate the
  /// faces
  /// @param color the (optional) color of the object to be written
  static void drawSurface(
      IVisualization& helper, const Surface& surface,
      const GeometryContext& gctx,
      const Transform3D& transform = Transform3D::Identity(), size_t lseg = 72,
      bool triangulate = false,
      const IVisualization::ColorType& color = {120, 120, 120});

  /// Helper method to Surface objects
  ///
  /// @param helper [in, out] The visualization helper
  /// @param surfaceArray The surface to be drawn
  /// @param gctx The geometry context for which it is drawn
  /// @param binning The binning prescritotion, if empty, the grid is not drawn
  /// @param transform An option additional transform
  /// @param lseg The number of segments for a full arch (if needed)
  /// @param triangulate The (optional) boolean diretive to triangulate the
  /// faces
  /// @param sfcolor the (optional) color for the surfaces
  /// @param gcolor the (otional) color of the grid
  static void drawSurfaceArray(
      IVisualization& helper, const SurfaceArray& surfaceArray,
      const GeometryContext& gctx, std::vector<BinningValue> binning = {},
      const Transform3D& transform = Transform3D::Identity(), size_t lseg = 72,
      bool triangulate = false,
      const IVisualization::ColorType& sfcolor = {120, 120, 120},
      const IVisualization::ColorType& gcolor = {200, 0, 0});

  /// Helper method for volume objects
  ///
  /// @param helper [in, out] The visualization helper
  /// @param volume The surface to be drawn
  /// @param gctx The geometry context for which it is drawn
  /// @param transform An option additional transform
  /// @param lseg The number of segments for a full arch (if needed)
  /// @param triangulate The (optional) boolean diretive to triangulate the
  /// faces
  /// @param color the (optional) color of the object to be written
  static void drawVolume(IVisualization& helper, const AbstractVolume& volume,
                         const GeometryContext& gctx,
                         const Transform3D& transform = Transform3D::Identity(),
                         size_t lseg = 72, bool triangulate = false,
                         const IVisualization::ColorType& color = {120, 120,
                                                                   120});

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
  static void drawSegmentBase(IVisualization& helper, const Vector3D& start,
                              const Vector3D& end, double thickness,
                              int arrows = 0, double arrowLength = 0.,
                              double arrowWidth = 0., size_t lseg = 72,
                              const IVisualization::ColorType& color = {
                                  120, 120, 120});

  /// Convenience function : line
  ///
  /// @param helper [in, out] The visualization helper
  /// @param start The start point
  /// @param end The end point
  /// @param thickness of the line, if bigger 0, approximated by cylinder
  /// @param lseg The number of segments for a full arch (if needed)
  /// @param color the (optional) color of the object to be written
  static void drawSegment(IVisualization& helper, const Vector3D& start,
                          const Vector3D& end, double thickness,
                          size_t lseg = 72,
                          const IVisualization::ColorType& color = {20, 120,
                                                                    120});

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
  static void drawArrowBackward(
      IVisualization& helper, const Vector3D& start, const Vector3D& end,
      double thickness, double arrowLength, double arrowWidth, size_t lseg = 72,
      const IVisualization::ColorType& color = {120, 120, 120});

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
  static void drawArrowForward(
      IVisualization& helper, const Vector3D& start, const Vector3D& end,
      double thickness, double arrowLength, double arrowWidth, size_t lseg = 72,
      const IVisualization::ColorType& color = {120, 120, 120});

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
  static void drawArrowsBoth(
      IVisualization& helper, const Vector3D& start, const Vector3D& end,
      double thickness, double arrowLength, double arrowWidth, size_t lseg = 72,
      const IVisualization::ColorType& color = {120, 120, 120});
};

}  // namespace Acts
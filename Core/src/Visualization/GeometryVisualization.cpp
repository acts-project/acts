// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Visualization/GeometryVisualization.hpp"
#include "Acts/Geometry/CylinderVolumeBounds.hpp"
#include "Acts/Geometry/Polyhedron.hpp"
#include "Acts/Surfaces/ConeBounds.hpp"
#include "Acts/Surfaces/ConeSurface.hpp"
#include "Acts/Surfaces/CylinderBounds.hpp"
#include "Acts/Surfaces/CylinderSurface.hpp"
#include "Acts/Surfaces/DiscSurface.hpp"
#include "Acts/Surfaces/RadialBounds.hpp"
#include "Acts/Utilities/UnitVectors.hpp"

void Acts::GeometryVisualization::drawSurface(
    IVisualization& helper, const Surface& surface, const GeometryContext& gctx,
    const Transform3D& transform, size_t lseg, bool triangulate,
    const IVisualization::ColorType& color) {
  // Drawing the polyhedron representation of surfaces
  Polyhedron phedron = surface.polyhedronRepresentation(gctx, lseg);
  phedron.move(transform);
  phedron.draw(helper, triangulate, color);
}

void Acts::GeometryVisualization::drawSurfaceArray(
    IVisualization& helper, const SurfaceArray& surfaceArray,
    const GeometryContext& gctx, std::vector<BinningValue> binning,
    const Transform3D& transform, size_t lseg, bool triangulate,
    const IVisualization::ColorType& sfcolor,
    const IVisualization::ColorType& gcolor) {
  // Draw all the surfaces
  Extent arrayExtent;
  for (const auto& sf : surfaceArray.surfaces()) {
    drawSurface(helper, *sf, gctx, transform, lseg, triangulate, sfcolor);
    auto sfExtent = sf->polyhedronRepresentation(gctx, 1).extent();
    arrayExtent.extend(sfExtent);
  }

  double thickness = 0.25;
  // Draw the grid itself
  auto axes = surfaceArray.getAxes();
  if (not binning.empty() and binning.size() == 2 and axes.size() == 2) {
    // Cylinder surface array
    if (binning[0] == binPhi and binning[1] == binZ) {
      double R = arrayExtent.medium(binR);
      auto phiValues = axes[0]->getBinEdges();
      auto zValues = axes[1]->getBinEdges();
      // Longitudinal lines
      for (auto phi : phiValues) {
        for (size_t iz = 1; iz < zValues.size(); ++iz) {
          double cphi = std::cos(phi);
          double sphi = std::sin(phi);
          Vector3D p1(R * cphi, R * sphi, zValues[iz]);
          Vector3D p0(R * cphi, R * sphi, zValues[iz - 1]);
          drawSegment(helper, p0, p1, 0.5 * thickness, 4, gcolor);
        }
        for (auto z : zValues) {
          CylinderVolumeBounds cvb(R - 0.5 * thickness, R + 0.5 * thickness,
                                   0.5 * thickness);
          auto cvbOrientedSurfaces = cvb.orientedSurfaces();
          for (auto cvbSf : cvbOrientedSurfaces) {
            drawSurface(helper, *cvbSf.first, gctx,
                        Translation3D(0., 0., z) * Transform3D::Identity(), 72,
                        triangulate, gcolor);
          }
        }
      }
    } else if (binning[0] == binR and binning[1] == binPhi) {
      double z = arrayExtent.medium(binZ);
      auto rValues = axes[0]->getBinEdges();
      auto phiValues = axes[1]->getBinEdges();
      for (auto r : rValues) {
        CylinderVolumeBounds cvb(r - 0.5 * thickness, r + 0.5 * thickness,
                                 0.5 * thickness);
        auto cvbOrientedSurfaces = cvb.orientedSurfaces();
        for (auto cvbSf : cvbOrientedSurfaces) {
          drawSurface(helper, *cvbSf.first, gctx,
                      Translation3D(0., 0., z) * Transform3D::Identity(), 72,
                      triangulate, gcolor);
        }
      }
      double rMin = axes[0]->getMin();
      double rMax = axes[0]->getMax();
      for (auto phi : phiValues) {
        double cphi = std::cos(phi);
        double sphi = std::sin(phi);
        Vector3D p1(rMax * cphi, rMax * sphi, z);
        Vector3D p0(rMin * cphi, rMin * sphi, z);
        drawSegment(helper, p0, p1, 0.5 * thickness, 4, gcolor);
      }
    }
  }
}

void Acts::GeometryVisualization::drawVolume(
    IVisualization& helper, const AbstractVolume& volume,
    const GeometryContext& gctx, const Transform3D& transform, size_t lseg,
    bool triangulate, const IVisualization::ColorType& color) {
  // Drawing the polyhedron representation of surfaces
  auto bSurfaces = volume.boundarySurfaces();
  for (const auto& bs : bSurfaces) {
    drawSurface(helper, bs->surfaceRepresentation(), gctx, transform, lseg,
                triangulate, color);
  }
}

void Acts::GeometryVisualization::drawSegmentBase(
    IVisualization& helper, const Vector3D& start, const Vector3D& end,
    double thickness, int arrows, double arrowLength, double arrowWidth,
    size_t lseg, const IVisualization::ColorType& color) {
  // Draw the parameter shaft and cone
  auto direction = Vector3D(end - start).normalized();
  double hlength = 0.5 * Vector3D(end - start).norm();

  auto unitVectors = makeCurvilinearUnitVectors(direction);
  RotationMatrix3D lrotation;
  lrotation.col(0) = unitVectors.first;
  lrotation.col(1) = unitVectors.second;
  lrotation.col(2) = direction;

  Vector3D lcenter = 0.5 * (start + end);
  double alength = (thickness > 0.) ? arrowLength * thickness : 2.;

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

void Acts::GeometryVisualization::drawSegment(
    IVisualization& helper, const Vector3D& start, const Vector3D& end,
    double thickness, size_t lseg, const IVisualization::ColorType& color) {
  drawSegmentBase(helper, start, end, thickness, 0, 0., 0., lseg, color);
}

void Acts::GeometryVisualization::drawArrowBackward(
    IVisualization& helper, const Vector3D& start, const Vector3D& end,
    double thickness, double arrowLength, double arrowWidth, size_t lseg,
    const IVisualization::ColorType& color) {
  drawSegmentBase(helper, start, end, thickness, -1, arrowLength, arrowWidth,
                  lseg, color);
}

void Acts::GeometryVisualization::drawArrowForward(
    IVisualization& helper, const Vector3D& start, const Vector3D& end,
    double thickness, double arrowLength, double arrowWidth, size_t lseg,
    const IVisualization::ColorType& color) {
  drawSegmentBase(helper, start, end, thickness, 1, arrowLength, arrowWidth,
                  lseg, color);
}

void Acts::GeometryVisualization::drawArrowsBoth(
    IVisualization& helper, const Vector3D& start, const Vector3D& end,
    double thickness, double arrowLength, double arrowWidth, size_t lseg,
    const IVisualization::ColorType& color) {
  drawSegmentBase(helper, start, end, thickness, 2, arrowLength, arrowWidth,
                  lseg, color);
}

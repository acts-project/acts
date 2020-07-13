// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Visualization/GeometryView.hpp"
#include "Acts/Geometry/CylinderVolumeBounds.hpp"
#include "Acts/Geometry/Layer.hpp"
#include "Acts/Geometry/Polyhedron.hpp"
#include "Acts/Geometry/TrackingVolume.hpp"
#include "Acts/Surfaces/ConeBounds.hpp"
#include "Acts/Surfaces/ConeSurface.hpp"
#include "Acts/Surfaces/CylinderBounds.hpp"
#include "Acts/Surfaces/CylinderSurface.hpp"
#include "Acts/Surfaces/DiscSurface.hpp"
#include "Acts/Surfaces/RadialBounds.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Surfaces/SurfaceArray.hpp"
#include "Acts/Utilities/UnitVectors.hpp"

void Acts::GeometryView::drawPolyhedron(IVisualization& helper,
                                        const Polyhedron& polyhedron,
                                        const ViewConfig& ViewConfig) {
  if (ViewConfig.visible) {
    if (not ViewConfig.triangulate) {
      helper.faces(polyhedron.vertices, polyhedron.faces, ViewConfig.color);
    } else {
      helper.faces(polyhedron.vertices, polyhedron.triangularMesh,
                   ViewConfig.color);
    }
  }
}

void Acts::GeometryView::drawSurface(IVisualization& helper,
                                     const Surface& surface,
                                     const GeometryContext& gctx,
                                     const Transform3D& transform,
                                     const ViewConfig& ViewConfig) {
  Polyhedron surfaceHedron =
      surface.polyhedronRepresentation(gctx, ViewConfig.nSegments);
  if (not transform.isApprox(Transform3D::Identity())) {
    surfaceHedron.move(transform);
  }
  drawPolyhedron(helper, surfaceHedron, ViewConfig);
}

void Acts::GeometryView::drawSurfaceArray(IVisualization& helper,
                                          const SurfaceArray& surfaceArray,
                                          const GeometryContext& gctx,
                                          const Transform3D& transform,
                                          const ViewConfig& sensitiveConfig,
                                          const ViewConfig& passiveConfig,
                                          const ViewConfig& gridConfig) {
  // Draw all the surfaces
  Extent arrayExtent;
  for (const auto& sf : surfaceArray.surfaces()) {
    ViewConfig vConfig = sf->associatedDetectorElement() != nullptr
                             ? sensitiveConfig
                             : passiveConfig;
    drawSurface(helper, *sf, gctx, transform, vConfig);
    auto sfExtent = sf->polyhedronRepresentation(gctx, 1).extent();
    arrayExtent.extend(sfExtent);
  }

  if (not sensitiveConfig.outputName.empty()) {
    helper.write(sensitiveConfig.outputName);
    helper.clear();
  }

  double thickness = gridConfig.lineThickness;
  // Draw the grid itself
  auto binning = surfaceArray.binningValues();
  auto axes = surfaceArray.getAxes();
  if (not binning.empty() and binning.size() == 2 and axes.size() == 2) {
    // Cylinder surface array
    if (binning[0] == binPhi and binning[1] == binZ) {
      double R = arrayExtent.medium(binR) + gridConfig.offset;
      auto phiValues = axes[0]->getBinEdges();
      auto zValues = axes[1]->getBinEdges();
      ViewConfig gridRadConfig = gridConfig;
      gridRadConfig.nSegments = phiValues.size();
      // Longitudinal lines
      for (auto phi : phiValues) {
        double cphi = std::cos(phi);
        double sphi = std::sin(phi);
        Vector3D p1(R * cphi, R * sphi, axes[1]->getMin());
        Vector3D p0(R * cphi, R * sphi, axes[1]->getMax());
        drawSegment(helper, transform * p0, transform * p1, gridConfig);
      }
      CylinderVolumeBounds cvb(R - 0.5 * thickness, R + 0.5 * thickness,
                               0.5 * thickness);
      auto cvbOrientedSurfaces = cvb.orientedSurfaces();
      for (auto z : zValues) {
        for (auto cvbSf : cvbOrientedSurfaces) {
          drawSurface(helper, *cvbSf.first, gctx,
                      Translation3D(0., 0., z) * transform, gridRadConfig);
        }
      }

    } else if (binning[0] == binR and binning[1] == binPhi) {
      double z = arrayExtent.medium(binZ) + gridConfig.offset;
      auto rValues = axes[0]->getBinEdges();
      auto phiValues = axes[1]->getBinEdges();
      ViewConfig gridRadConfig = gridConfig;
      gridRadConfig.nSegments = phiValues.size();
      for (auto r : rValues) {
        CylinderVolumeBounds cvb(r - 0.5 * thickness, r + 0.5 * thickness,
                                 0.5 * thickness);
        auto cvbOrientedSurfaces = cvb.orientedSurfaces();
        for (auto cvbSf : cvbOrientedSurfaces) {
          drawSurface(helper, *cvbSf.first, gctx,
                      Translation3D(0., 0., z) * transform, gridRadConfig);
        }
      }
      double rMin = axes[0]->getMin();
      double rMax = axes[0]->getMax();
      for (auto phi : phiValues) {
        double cphi = std::cos(phi);
        double sphi = std::sin(phi);
        Vector3D p1(rMax * cphi, rMax * sphi, z);
        Vector3D p0(rMin * cphi, rMin * sphi, z);
        drawSegment(helper, transform * p0, transform * p1, gridConfig);
      }
    }
  }

  if (not gridConfig.outputName.empty()) {
    helper.write(gridConfig.outputName);
    helper.clear();
  }
}

void Acts::GeometryView::drawVolume(IVisualization& helper,
                                    const AbstractVolume& volume,
                                    const GeometryContext& gctx,
                                    const Transform3D& transform,
                                    const ViewConfig& viewConfig) {
  auto bSurfaces = volume.boundarySurfaces();
  for (const auto& bs : bSurfaces) {
    drawSurface(helper, bs->surfaceRepresentation(), gctx, transform,
                viewConfig);
  }
}

void Acts::GeometryView::drawLayer(IVisualization& helper, const Layer& layer,
                                   const GeometryContext& gctx,
                                   const ViewConfig& layerConfig,
                                   const ViewConfig& sensitiveConfig,
                                   const ViewConfig& gridConfig) {
  if (layerConfig.visible) {
    auto layerVolume = layer.representingVolume();
    if (layerVolume != nullptr) {
      drawVolume(helper, *layerVolume, gctx, Transform3D::Identity(),
                 layerConfig);
    } else {
      const auto& layerSurface = layer.surfaceRepresentation();
      drawSurface(helper, layerSurface, gctx, Transform3D::Identity(),
                  layerConfig);
    }
    if (not layerConfig.outputName.empty()) {
      helper.write(layerConfig.outputName);
      helper.clear();
    }
  }

  if (sensitiveConfig.visible or gridConfig.visible) {
    auto surfaceArray = layer.surfaceArray();
    if (surfaceArray != nullptr) {
      drawSurfaceArray(helper, *surfaceArray, gctx, Transform3D::Identity(),
                       sensitiveConfig, layerConfig, gridConfig);
    }
  }
}

void Acts::GeometryView::drawTrackingVolume(
    IVisualization& helper, const TrackingVolume& tVolume,
    const GeometryContext& gctx, const ViewConfig& containerView,
    const ViewConfig& volumeView, const ViewConfig& layerView,
    const ViewConfig& sensitiveView, const ViewConfig& gridView, bool writeIt,
    const std::string& tag) {
  if (tVolume.confinedVolumes() != nullptr) {
    const auto& subVolumes = tVolume.confinedVolumes()->arrayObjects();
    for (const auto& tv : subVolumes) {
      drawTrackingVolume(helper, *tv, gctx, containerView, volumeView,
                         layerView, sensitiveView, gridView, writeIt, tag);
    }
  }

  ViewConfig cConfig = containerView;
  ViewConfig vConfig = volumeView;
  ViewConfig lConfig = layerView;
  ViewConfig sConfig = sensitiveView;
  ViewConfig gConfig = gridView;
  gConfig.nSegments = 8;

  ViewConfig vcConfig = cConfig;
  std::string vname = tVolume.volumeName();
  if (writeIt) {
    std::vector<std::string> repChar = {":", "|", " ", "{", "}"};
    for (auto rchar : repChar) {
      while (vname.find(rchar) != std::string::npos) {
        vname.replace(vname.find(rchar), rchar.size(), std::string("_"));
      }
    }
    if (tVolume.confinedVolumes() == nullptr) {
      vcConfig = vConfig;
      vcConfig.outputName = vname + std::string("_boundaries") + tag;
    } else {
      vcConfig.outputName =
          std::string("Container-") + vname + std::string("_boundaries") + tag;
    }
  }

  auto bSurfaces = tVolume.boundarySurfaces();
  for (const auto& bs : bSurfaces) {
    drawSurface(helper, bs->surfaceRepresentation(), gctx,
                Transform3D::Identity(), vcConfig);
  }
  if (writeIt) {
    helper.write(vcConfig.outputName);
    helper.clear();
  }

  if (tVolume.confinedLayers() != nullptr) {
    const auto& layers = tVolume.confinedLayers()->arrayObjects();
    size_t il = 0;
    for (const auto& tl : layers) {
      if (writeIt) {
        lConfig.outputName =
            vname + std::string("_passives_l") + std::to_string(il) + tag;
        sConfig.outputName =
            vname + std::string("_sensitives_l") + std::to_string(il) + tag;
        gConfig.outputName =
            vname + std::string("_grids_l") + std::to_string(il) + tag;
      }
      drawLayer(helper, *tl, gctx, lConfig, sConfig, gConfig);
      ++il;
    }
  }
}

void Acts::GeometryView::drawSegmentBase(IVisualization& helper,
                                         const Vector3D& start,
                                         const Vector3D& end, int arrows,
                                         double arrowLength, double arrowWidth,
                                         const ViewConfig& viewConfig) {
  double thickness = viewConfig.lineThickness;

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
  if (alength > hlength) {
    alength = hlength;
  }

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

    drawSurface(helper, *line, GeometryContext(), Transform3D::Identity(),
                viewConfig);
  } else {
    helper.line(start, end, viewConfig.color);
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
                  viewConfig);
      // Arrow end plate
      auto aptransform = std::make_shared<Transform3D>(Transform3D::Identity());
      aptransform->prerotate(lrotation);
      aptransform->pretranslate(Vector3D(end - alength * direction));

      auto plate = Surface::makeShared<DiscSurface>(aptransform, plateBounds);
      drawSurface(helper, *plate, GeometryContext(), Transform3D::Identity(),
                  viewConfig);
    }
    if (arrows < 0 or arrows == 2) {
      auto astransform = std::make_shared<Transform3D>(Transform3D::Identity());
      astransform->prerotate(lrotation);
      astransform->pretranslate(start);

      // Arrow cone
      auto coneBounds = std::make_shared<ConeBounds>(alpha, 0., alength);
      auto cone = Surface::makeShared<ConeSurface>(astransform, coneBounds);
      drawSurface(helper, *cone, GeometryContext(), Transform3D::Identity(),
                  viewConfig);
      // Arrow end plate
      auto aptransform = std::make_shared<Transform3D>(Transform3D::Identity());
      aptransform->prerotate(lrotation);
      aptransform->pretranslate(Vector3D(start + alength * direction));

      auto plate = Surface::makeShared<DiscSurface>(aptransform, plateBounds);
      drawSurface(helper, *plate, GeometryContext(), Transform3D::Identity(),
                  viewConfig);
    }
  }
}

void Acts::GeometryView::drawSegment(IVisualization& helper,
                                     const Vector3D& start, const Vector3D& end,
                                     const ViewConfig& viewConfig) {
  drawSegmentBase(helper, start, end, 0, 0., 0., viewConfig);
}

void Acts::GeometryView::drawArrowBackward(
    IVisualization& helper, const Vector3D& start, const Vector3D& end,
    double arrowLength, double arrowWidth, const ViewConfig& viewConfig) {
  drawSegmentBase(helper, start, end, -1, arrowLength, arrowWidth, viewConfig);
}

void Acts::GeometryView::drawArrowForward(IVisualization& helper,
                                          const Vector3D& start,
                                          const Vector3D& end,
                                          double arrowLength, double arrowWidth,
                                          const ViewConfig& viewConfig) {
  drawSegmentBase(helper, start, end, 1, arrowLength, arrowWidth, viewConfig);
}

void Acts::GeometryView::drawArrowsBoth(IVisualization& helper,
                                        const Vector3D& start,
                                        const Vector3D& end, double arrowLength,
                                        double arrowWidth,
                                        const ViewConfig& viewConfig) {
  drawSegmentBase(helper, start, end, 2, arrowLength, arrowWidth, viewConfig);
}

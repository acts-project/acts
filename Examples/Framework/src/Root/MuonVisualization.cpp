// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/Root/MuonVisualization.hpp"

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/Units.hpp"
#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "Acts/Geometry/TrackingVolume.hpp"
#include "Acts/Surfaces/LineBounds.hpp"
#include "Acts/Surfaces/RectangleBounds.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Surfaces/SurfaceBounds.hpp"
#include "Acts/Surfaces/TrapezoidBounds.hpp"
#include "Acts/Utilities/ArrayHelpers.hpp"
#include "Acts/Utilities/RangeXD.hpp"
#include "ActsExamples/EventData/SimParticle.hpp"

#include <cassert>
#include <memory>
#include <vector>

#include "TArrow.h"
#include "TBox.h"
#include "TCanvas.h"
#include "TEllipse.h"
#include "TH2I.h"
#include "TObject.h"

using namespace Acts;
using namespace Acts::UnitLiterals;

namespace {

/// @brief Returns the chamber ID from a full geometry identifier
GeometryIdentifier toChamberId(const GeometryIdentifier& id) {
  return GeometryIdentifier{}.withVolume(id.volume()).withLayer(id.layer());
}

/// @brief Returns the half-height of a trapezoid / rectangular bounds
double halfHeight(const SurfaceBounds& bounds) {
  if (bounds.type() == SurfaceBounds::BoundsType::eRectangle) {
    return static_cast<const RectangleBounds&>(bounds).get(
        RectangleBounds::eMaxY);
  }
  // Trapezoid -> endcap
  else if (bounds.type() == SurfaceBounds::BoundsType::eTrapezoid) {
    return static_cast<const TrapezoidBounds&>(bounds).get(
        TrapezoidBounds::eHalfLengthY);
  }
  return std::numeric_limits<double>::max();
}

/// @brief Draw a circle at position
std::unique_ptr<TEllipse> drawCircle(const double x, double y, const double r,
                                     const int color = kBlack,
                                     const int fillStyle = 0) {
  auto circle = std::make_unique<TEllipse>(x, y, r);
  circle->SetLineColor(color);
  circle->SetFillStyle(fillStyle);
  circle->SetLineWidth(1);
  circle->SetFillColorAlpha(color, 0.2);
  return circle;
}

/// @brief Draw a box at position
std::unique_ptr<TBox> drawBox(const double cX, const double cY, const double wX,
                              const double wY, const int color = kBlack,
                              const int fillStyle = 3344) {
  auto box = std::make_unique<TBox>(cX - 0.5 * wX, cY - 0.5 * wY, cX + 0.5 * wX,
                                    cY + 0.5 * wY);
  box->SetLineColor(color);
  box->SetFillStyle(fillStyle);
  box->SetLineWidth(1);
  box->SetFillColorAlpha(color, 0.2);
  return box;
}

/// @brief Determines the axis ranges in the z-y plane to draw
///        the measurements on the canvas
RangeXD<2, double> canvasRanges(
    const ActsExamples::MuonSpacePointBucket& bucket) {
  RangeXD<2, double> ranges{
      filledArray<double, 2>(std::numeric_limits<double>::max()),
      filledArray<double, 2>(-std::numeric_limits<double>::max())};
  // Extra margin such that the canvas axes don't overlap with the depicted
  // measurements
  constexpr double extra = 3._cm;
  for (const auto& sp : bucket) {
    const Vector3& pos = sp.localPosition();
    ranges.expand(0, pos.z() - extra, pos.z() + extra);
    ranges.expand(1, pos.y() - extra, pos.y() + extra);
  }
  return ranges;
}

/// @brief Returns whether a surface should be drawn on the canvas
bool isSurfaceToDraw(
    const Acts::GeometryContext& gctx, const Surface& surface,
    const RangeXD<2, double>& canvasBoundaries,
    const std::function<Acts::Transform3(const Acts::GeometryContext&,
                                         const Acts::GeometryIdentifier&)>&
        toSpacePointFrame) {
  // Draw only active surfaces
  if (surface.associatedDetectorElement() == nullptr) {
    return false;
  }
  // surface position in the frame
  const Vector3 pos =
      toSpacePointFrame(gctx, surface.geometryId()).translation();
  const auto& bounds = surface.bounds();

  if (surface.type() == Surface::SurfaceType::Plane) {
    const double hL{halfHeight(bounds)};
    // check whether the surface is inside the visible range
    const double minZ = std::max(pos.z() - hL, canvasBoundaries.min(0));
    const double maxZ = std::min(pos.z() + hL, canvasBoundaries.max(0));
    // The maximum is below the left side of the strip plane
    // or the minimum is above the right side of the plane
    return maxZ > pos.z() - hL && minZ < pos.z() + hL &&
           pos.y() > canvasBoundaries.min(1) &&
           pos.y() < canvasBoundaries.max(1);
  } else if (surface.type() == Surface::SurfaceType::Straw) {
    const double r = static_cast<const LineBounds&>(bounds).get(LineBounds::eR);
    // Check that the straw surface is well embedded on the canvas
    return pos.y() - r > canvasBoundaries.min(1) &&
           pos.y() + r < canvasBoundaries.max(1) &&
           pos.z() - r > canvasBoundaries.min(0) &&
           pos.z() + r < canvasBoundaries.max(0);
  }

  return false;
}

}  // namespace

namespace ActsExamples {

void visualizeMuonSpacePoints(
    const std::string& outputPath, const GeometryContext& gctx,
    const MuonSpacePointBucket& bucket, const SimHitContainer& simHits,
    const SimParticleContainer& simParticles,
    const std::function<Transform3(
        const GeometryContext&, const GeometryIdentifier&)>& toSpacePointFrame,
    const TrackingGeometry& trackingGeometry) {
  auto canvas = std::make_unique<TCanvas>("can", "can", 600, 600);
  canvas->cd();

  const GeometryIdentifier chambId = toChamberId(bucket.front().geometryId());

  std::vector<std::unique_ptr<TObject>> primitives{};

  const RangeXD<2, double> canvasBound{canvasRanges(bucket)};
  /// Draw the frame
  auto frame = std::make_unique<TH2I>("frame", "frame;z [mm];y [mm]", 1,
                                      canvasBound.min(0), canvasBound.max(0), 1,
                                      canvasBound.min(1), canvasBound.max(1));
  frame->Draw("AXIS");

  // Loop over all surfaces inside the chamber volume to draw the ones covered
  // by the canvas
  const TrackingVolume* chambVolume = trackingGeometry.findVolume(chambId);
  assert(chambVolume != nullptr);
  chambVolume->apply(overloaded{[&canvasBound, &gctx, &primitives,
                                 &toSpacePointFrame](const Surface& surface) {
    if (!isSurfaceToDraw(gctx, surface, canvasBound, toSpacePointFrame)) {
      return;
    }
    const Vector3 pos =
        toSpacePointFrame(gctx, surface.geometryId()).translation();
    const auto& bounds = surface.bounds();
    if (surface.type() == Surface::SurfaceType::Plane) {
      const double hL{halfHeight(bounds)};
      const double minZ = std::max(pos.z() - hL, canvasBound.min(0));
      const double maxZ = std::min(pos.z() + hL, canvasBound.max(0));
      primitives.push_back(drawBox(0.5 * (minZ + maxZ), pos.y(), maxZ - minZ,
                                   0.3_cm, kBlack, 0));

    } else if (surface.type() == Surface::SurfaceType::Straw) {
      const double r =
          static_cast<const LineBounds&>(bounds).get(LineBounds::eR);
      primitives.push_back(drawCircle(pos.z(), pos.y(), r, kBlack, 0));
    }
  }});

  for (auto& sp : bucket) {
    const Vector3& pos = sp.localPosition();
    if (sp.isStraw()) {
      primitives.push_back(
          drawCircle(pos.z(), pos.y(), sp.driftRadius(), kRed, 0));
    } else {
      primitives.push_back(
          drawBox(pos.z(), pos.y(), 3._cm, 0.5_cm, kRed, 1001));
    }
  }
  // Finally draw the muon trajectory
  for (const auto& simHit : simHits) {
    if (chambId != toChamberId(simHit.geometryId())) {
      continue;
    }
    const auto simPartItr = simParticles.find(simHit.particleId());
    if (simPartItr == simParticles.end() ||
        (*simPartItr).hypothesis() != ParticleHypothesis::muon()) {
      continue;
    }
    const auto toSpTrf = toSpacePointFrame(gctx, simHit.geometryId()) *
                         trackingGeometry.findSurface(simHit.geometryId())
                             ->localToGlobalTransform(gctx)
                             .inverse();
    const Vector3 pos = toSpTrf * simHit.position();
    const Vector3 dir = toSpTrf.linear() * simHit.direction();
    constexpr double arrowL = 1._cm;
    const Vector3 start = pos - 0.5 * arrowL * dir;
    const Vector3 end = pos + 0.5 * arrowL * dir;
    auto arrow =
        std::make_unique<TArrow>(start.z(), start.y(), end.z(), end.y(), 0.03);
    arrow->SetLineColor(kBlue + 1);
    arrow->SetLineWidth(1);
    arrow->SetLineStyle(kSolid);
    primitives.push_back(std::move(arrow));
  }

  for (auto& prim : primitives) {
    prim->Draw();
  }
  canvas->SaveAs(outputPath.c_str());
}

}  // namespace ActsExamples

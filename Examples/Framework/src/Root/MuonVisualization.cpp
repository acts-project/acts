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
#include "Acts/Utilities/Logger.hpp"
#include "Acts/Utilities/RangeXD.hpp"
#include "Acts/Utilities/StringHelpers.hpp"
#include "ActsExamples/EventData/SimParticle.hpp"

#include <algorithm>
#include <cassert>
#include <format>
#include <memory>
#include <mutex>
#include <vector>

#include "TArrow.h"
#include "TBox.h"
#include "TCanvas.h"
#include "TEllipse.h"
#include "TH2D.h"
#include "TH2I.h"
#include "TLegend.h"
#include "TMarker.h"
#include "TObject.h"
#include "TStyle.h"

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
bool isSurfaceToDraw(const Surface& surface,
                     const RangeXD<2, double>& canvasBoundaries,
                     const Vector3& posInFrame) {
  // Draw only active surfaces
  if (surface.surfacePlacement() == nullptr) {
    return false;
  }
  const auto& bounds = surface.bounds();

  if (surface.type() == Surface::SurfaceType::Plane) {
    const double hL{halfHeight(bounds)};
    // check whether the surface is inside the visible range
    const double minZ = std::max(posInFrame.z() - hL, canvasBoundaries.min(0));
    const double maxZ = std::min(posInFrame.z() + hL, canvasBoundaries.max(0));
    // The maximum is below the left side of the strip plane
    // or the minimum is above the right side of the plane
    return maxZ > posInFrame.z() - hL && minZ < posInFrame.z() + hL &&
           posInFrame.y() > canvasBoundaries.min(1) &&
           posInFrame.y() < canvasBoundaries.max(1);
  } else if (surface.type() == Surface::SurfaceType::Straw) {
    const double r = static_cast<const LineBounds&>(bounds).get(LineBounds::eR);
    // Check that the straw surface is well embedded on the canvas
    return posInFrame.y() - r > canvasBoundaries.min(1) &&
           posInFrame.y() + r < canvasBoundaries.max(1) &&
           posInFrame.z() - r > canvasBoundaries.min(0) &&
           posInFrame.z() + r < canvasBoundaries.max(0);
  }

  return false;
}

}  // namespace

namespace ActsExamples {

void visualizeMuonSpacePoints(const std::string& outputPath,
                              const GeometryContext& gctx,
                              const MuonSpacePointBucket& bucket,
                              const SimHitContainer& simHits,
                              const SimParticleContainer& simParticles,
                              const TrackingGeometry& trackingGeometry,
                              const Acts::Logger& logger) {
  // Helper function to transform from hit frame to spacepoint frame
  auto toSpacePointFrame = [&](const GeometryIdentifier& hitId) -> Transform3 {
    const Surface* hitSurf = trackingGeometry.findSurface(hitId);
    assert(hitSurf != nullptr);

    // Fetch the parent volume to express all points in the same coordinate
    // system
    const TrackingVolume* volume =
        trackingGeometry.findVolume(toChamberId(hitId));
    assert(volume != nullptr);

    // Transformation to the common coordinate system of all spacepoints
    const Transform3 parentTrf{AngleAxis3{90._degree, Vector3::UnitZ()} *
                               volume->globalToLocalTransform(gctx) *
                               hitSurf->localToGlobalTransform(gctx)};
    ACTS_VERBOSE("Transform into spacepoint frame for surface "
                 << hitId << " is \n"
                 << Acts::toString(parentTrf));
    return parentTrf;
  };

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
  chambVolume->apply(overloaded{[&](const Surface& surface) {
    const Vector3 pos = toSpacePointFrame(surface.geometryId()).translation();
    if (!isSurfaceToDraw(surface, canvasBound, pos)) {
      return;
    }
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
    const auto toSpTrf = toSpacePointFrame(simHit.geometryId()) *
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

void visualizeMuonHoughMaxima(
    const std::string& outputPath, const MuonSpacePoint::MuonId& bucketId,
    const std::vector<Acts::HoughTransformUtils::PeakFinders::IslandsAroundMax<
        const MuonSpacePoint*>::Maximum>& maxima,
    const Acts::HoughTransformUtils::HoughPlane<const MuonSpacePoint*>& plane,
    const Acts::HoughTransformUtils::HoughAxisRanges& axis,
    const MuonSegmentContainer& truthSegments, const Acts::Logger& logger) {
  static std::mutex canvasMutex{};
  std::lock_guard guard{canvasMutex};

  TCanvas canvas("houghCanvas", "", 800, 800);
  canvas.SetRightMargin(0.12);
  canvas.SetLeftMargin(0.12);
  gStyle->SetPalette(kGreyScale);
  gStyle->SetOptStat(0);

  std::vector<std::unique_ptr<TObject>> primitives;

  /// Save the hough accumulator as histogram
  TH2D houghHistoForPlot("houghHist", "HoughPlane;tan(#alpha);z0 [mm]",
                         plane.nBinsX(), axis.xMin, axis.xMax, plane.nBinsY(),
                         axis.yMin, axis.yMax);
  houghHistoForPlot.SetTitle(
      std::format("Station {:}, side {:}, sector {:2d}",
                  MuonSpacePoint::MuonId::toString(bucketId.msStation()),
                  MuonSpacePoint::MuonId::toString(bucketId.side()),
                  bucketId.sector())
          .c_str());

  /** Copy the plane content into the histogram */
  for (int bx = 0; bx < houghHistoForPlot.GetNbinsX(); ++bx) {
    for (int by = 0; by < houghHistoForPlot.GetNbinsY(); ++by) {
      houghHistoForPlot.SetBinContent(bx + 1, by + 1, plane.nHits(bx, by));
    }
  }
  /** Set the contours */
  auto maxHitsAsInt = static_cast<int>(plane.maxHits());
  houghHistoForPlot.SetContour(maxHitsAsInt + 1);
  for (int k = 0; k < maxHitsAsInt + 1; ++k) {
    houghHistoForPlot.SetContourLevel(k, k - 0.5);
  }

  auto legend = std::make_unique<TLegend>(0.5, 0.7, 1. - gPad->GetRightMargin(),
                                          1. - gPad->GetTopMargin());
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);

  /** Fetch the true parameters */
  MuonSegmentContainer::const_iterator truthItr = truthSegments.begin();
  const int trueCoord = bucketId.measuresEta() ? Acts::eY : Acts::eX;
  while ((truthItr = std::find_if(truthItr, truthSegments.end(),
                                  [bucketId](const MuonSegment& seg) {
                                    return seg.id().sameStation(bucketId);
                                  })) != truthSegments.end()) {
    const MuonSegment& truthSeg{*truthItr};
    const float tanAlpha =
        truthSeg.localDirection()[trueCoord] / truthSeg.localDirection().z();
    const float intercept = truthSeg.localPosition()[trueCoord];

    auto trueMarker =
        std::make_unique<TMarker>(tanAlpha, intercept, kOpenCrossX);
    trueMarker->SetMarkerSize(3);
    trueMarker->SetMarkerColor(kRed);
    legend->AddEntry(trueMarker.get(), "True coordinates");
    primitives.push_back(std::move(trueMarker));
    ACTS_VERBOSE("Draw true segment " << truthSeg.id()
                                      << " in bucket: " << bucketId);
    ++truthItr;
  }

  bool addedLeg{false};
  for (const auto& max : maxima) {
    auto marker = std::make_unique<TMarker>(max.x, max.y, kFullSquare);
    marker->SetMarkerSize(1);
    marker->SetMarkerColor(kBlue);

    auto box = std::make_unique<TBox>(max.x - max.wx, max.y - max.wy,
                                      max.x + max.wx, max.y + max.wy);
    box->SetLineColor(kBlue);
    box->SetFillStyle(1001);
    box->SetFillColorAlpha(kBlue, 0.1);
    box->SetLineWidth(0);
    if (!addedLeg) {
      legend->AddEntry(marker.get(), "Hough maxima");
      legend->AddEntry(box.get(), "Hough uncertainties");
      addedLeg = true;
    }
    primitives.push_back(std::move(box));
    primitives.push_back(std::move(marker));
  }
  primitives.emplace_back(std::move(legend));
  for (auto& prim : primitives) {
    prim->Draw();
  }

  canvas.SaveAs(outputPath.c_str());
}

}  // namespace ActsExamples

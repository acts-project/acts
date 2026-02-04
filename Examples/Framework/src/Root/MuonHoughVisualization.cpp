// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/Root/MuonHoughVisualization.hpp"

#include "Acts/Definitions/Units.hpp"

#include <algorithm>
#include <format>
#include <memory>
#include <mutex>
#include <vector>

#include "TBox.h"
#include "TCanvas.h"
#include "TH2D.h"
#include "TLegend.h"
#include "TMarker.h"
#include "TStyle.h"

namespace ActsExamples {

void visualizeMuonHoughMaxima(
    const std::string& outputPath, const MuonSpacePoint::MuonId& bucketId,
    const std::vector<Acts::HoughTransformUtils::PeakFinders::
                          IslandsAroundMax<const MuonSpacePoint*>::Maximum>&
        maxima,
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

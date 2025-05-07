// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/TrackFinding/MuonHoughSeeder.hpp"

#include "Acts/Definitions/Units.hpp"
#include "ActsExamples/EventData/MuonHoughMaximum.hpp"
#include "ActsExamples/EventData/MuonSegment.hpp"
#include "ActsExamples/EventData/MuonSpacePoint.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

#include "TBox.h"
#include "TCanvas.h"
#include "TH2D.h"
#include "TLegend.h"
#include "TMarker.h"
#include "TStyle.h"

namespace ActsExamples {

// left solution
auto etaHoughParamDC_left = [](double tanAlpha, const MuonSpacePoint& DC) {
  return DC.localPosition().y() - tanAlpha * DC.localPosition().z() -
         DC.driftRadius() * std::sqrt(1. + tanAlpha * tanAlpha);
};
// right solution
auto etaHoughParamDC_right = [](double tanAlpha, const MuonSpacePoint& DC) {
  return DC.localPosition().y() - tanAlpha * DC.localPosition().z() +
         DC.driftRadius() * std::sqrt(1. + tanAlpha * tanAlpha);
};

// create the function parametrising the drift radius uncertainty
auto houghWidth_fromDC = [](double, const MuonSpacePoint& DC) {
  // scale reported errors up to at least 1mm or 3 times the reported error as
  // drift circle calib not fully reliable at this stage
  return std::min(std::sqrt(DC.covariance()(Acts::eY, Acts::eY)) * 3., 1.0);
};

/// strip solution
auto etaHoughParam_strip = [](double tanAlpha, const MuonSpacePoint& strip) {
  return strip.localPosition().y() - tanAlpha * strip.localPosition().z();
};
auto phiHoughParam_strip = [](double tanBeta, const MuonSpacePoint& strip) {
  return strip.localPosition().x() - tanBeta * strip.localPosition().z();
};

/// @brief Strip uncertainty
auto etaHoughWidth_strip = [](double, const MuonSpacePoint& strip) {
  return std::sqrt(strip.covariance()(Acts::eY, Acts::eY)) * 3.;
};
auto phiHoughWidth_strip = [](double, const MuonSpacePoint& strip) {
  return std::sqrt(strip.covariance()(Acts::eX, Acts::eX)) * 3.;
};

MuonHoughSeeder::MuonHoughSeeder(MuonHoughSeeder::Config cfg,
                                 Acts::Logging::Level lvl)
    : IAlgorithm("MuonHoughSeeder", lvl),
      m_cfg(std::move(cfg)),
      m_logger(Acts::getDefaultLogger(name(), lvl)) {
  if (m_cfg.inSpacePoints.empty()) {
    throw std::invalid_argument(
        "MuonHoughSeeder: Missing drift circle collection");
  }
  if (m_cfg.inTruthSegments.empty()) {
    throw std::invalid_argument(
        "MuonHoughSeeder: Missing truth segment collection");
  }

  m_inputSpacePoints.initialize(m_cfg.inSpacePoints);
  m_inputTruthSegs.initialize(m_cfg.inTruthSegments);
  m_outputMaxima.initialize(m_cfg.outHoughMax);
}

MuonHoughSeeder::~MuonHoughSeeder() = default;

ProcessCode MuonHoughSeeder::execute(const AlgorithmContext& ctx) const {
  // read the hits and circles
  const MuonSpacePointContainer& gotSpacePoints = m_inputSpacePoints(ctx);

  // configure the binning of the hough plane
  Acts::HoughTransformUtils::HoughPlaneConfig etaPlaneCfg;
  etaPlaneCfg.nBinsX = 25;
  etaPlaneCfg.nBinsY = 25;

  Acts::HoughTransformUtils::HoughPlaneConfig phiPlaneCfg;
  phiPlaneCfg.nBinsX = 10;
  phiPlaneCfg.nBinsY = 10;

  // instantiate the hough plane
  HoughPlane_t etaPlane{etaPlaneCfg};
  HoughPlane_t phiPlane{phiPlaneCfg};
  MuonHoughMaxContainer outMaxima{};

  for (const MuonSpacePointContainer::value_type& bucket : gotSpacePoints) {
    MuonHoughMaxContainer etaMaxima = constructEtaMaxima(ctx, bucket, etaPlane);
    ACTS_VERBOSE(__func__ << "() - " << __LINE__ << " Found "
                          << etaMaxima.size() << " eta maxima");
    MuonHoughMaxContainer stationMaxima =
        extendMaximaWithPhi(ctx, std::move(etaMaxima), phiPlane);
    ACTS_VERBOSE(__func__ << "() - " << __LINE__ << " Extended the maxima to "
                          << stationMaxima.size() << " station maxima");
    outMaxima.insert(outMaxima.end(),
                     std::make_move_iterator(stationMaxima.begin()),
                     std::make_move_iterator(stationMaxima.end()));
  }

  m_outputMaxima(ctx, std::move(outMaxima));

  return ProcessCode::SUCCESS;
}

MuonHoughMaxContainer MuonHoughSeeder::constructEtaMaxima(
    const AlgorithmContext& ctx, const MuonSpacePointBucket& bucket,
    HoughPlane_t& plane) const {
  MuonHoughMaxContainer etaMaxima{};
  AxisRange_t axisRanges{-3., 3., 100. * Acts::UnitConstants::m,
                         -100. * Acts::UnitConstants::m};
  /** brief Adapt the intercept axis ranges */
  for (const MuonSpacePoint& sp : bucket) {
    /** Require that the space point measures eta */
    if (!sp.id().measuresEta()) {
      continue;
    }
    const double y = sp.localPosition().y();

    axisRanges.yMin = std::min(axisRanges.yMin, y - m_cfg.etaPlaneMarginIcept);
    axisRanges.yMax = std::min(axisRanges.yMax, y + m_cfg.etaPlaneMarginIcept);
  }
  /** Ranges are adapted. Now fill the hough plane */
  plane.reset();
  for (const MuonSpacePoint& sp : bucket) {
    /** Require that the space point measures eta */
    if (!sp.id().measuresEta()) {
      continue;
    }
    if (sp.id().technology() == MuonSpacePoint::MuonId::TechField::Mdt) {
      plane.fill<MuonSpacePoint>(sp, axisRanges, etaHoughParamDC_left,
                                 houghWidth_fromDC, &sp, sp.id().detLayer());
      plane.fill<MuonSpacePoint>(sp, axisRanges, etaHoughParamDC_right,
                                 houghWidth_fromDC, &sp, sp.id().detLayer());
    } else {
      plane.fill<MuonSpacePoint>(sp, axisRanges, etaHoughParam_strip,
                                 etaHoughWidth_strip, &sp, sp.id().detLayer());
    }
  }
  /** Extract the maxima from the peak */
  PeakFinderCfg_t peakFinderCfg{};
  peakFinderCfg.fractionCutoff = 0.7f;
  peakFinderCfg.threshold = 3.f;
  peakFinderCfg.minSpacingBetweenPeaks = {0.f, 30.f};

  PeakFinder_t peakFinder{peakFinderCfg};

  /// Fetch the first set of peaks
  const MaximumVec_t maxima = peakFinder.findPeaks(plane, axisRanges);
  if (maxima.empty()) {
    ACTS_VERBOSE(__func__ << "() - " << __LINE__
                          << " No maxima found for bucket");
    return etaMaxima;
  }
  /** Create a measurement vector of the pure phi hits to be copied on all eta
   * maxima */
  MuonHoughMaximum::HitVec phiHits{};
  for (const MuonSpacePoint& sp : bucket) {
    if (!sp.id().measuresEta()) {
      phiHits.emplace_back(&sp);
    }
  }
  for (const Maximum_t& peakMax : maxima) {
    MuonHoughMaximum::HitVec hits{peakMax.hitIdentifiers.begin(),
                                  peakMax.hitIdentifiers.end()};
    hits.insert(hits.end(), phiHits.begin(), phiHits.end());
    etaMaxima.emplace_back(peakMax.x, peakMax.y, std::move(hits));
  }
  MuonId bucketId{bucket.front().id()};
  bucketId.setCoordFlags(true, false);
  displayMaxima(ctx, bucketId, maxima, plane, axisRanges);
  return etaMaxima;
}
MuonHoughMaxContainer MuonHoughSeeder::extendMaximaWithPhi(
    const AlgorithmContext& ctx, MuonHoughMaxContainer&& etaMaxima,
    HoughPlane_t& plane) const {
  MuonHoughMaxContainer outMaxima{};

  PeakFinderCfg_t peakFinderCfg{};
  peakFinderCfg.fractionCutoff = 0.7;
  peakFinderCfg.threshold = 2.;
  peakFinderCfg.minSpacingBetweenPeaks = {0., 30.};

  PeakFinder_t peakFinder{peakFinderCfg};

  for (MuonHoughMaximum& etaMax : etaMaxima) {
    plane.reset();
    AxisRange_t axisRanges{-3., 3., 100. * Acts::UnitConstants::m,
                           -100. * Acts::UnitConstants::m};

    unsigned nPhiHits{0};
    for (const MuonSpacePoint* sp : etaMax.hits()) {
      if (!sp->id().measuresPhi()) {
        continue;
      }
      ++nPhiHits;
      const double x = sp->localPosition().x();
      axisRanges.yMax =
          std::min(axisRanges.yMax, x + m_cfg.phiPlaneMarginIcept);
      axisRanges.yMin =
          std::min(axisRanges.yMin, x - m_cfg.phiPlaneMarginIcept);
    }
    if (nPhiHits < 2) {
      outMaxima.emplace_back(std::move(etaMax));
      continue;
    }
    MuonHoughMaximum::HitVec etaHits{};
    etaHits.reserve(etaMax.hits().size());
    for (const MuonSpacePoint* sp : etaMax.hits()) {
      if (!sp->id().measuresPhi()) {
        etaHits.push_back(sp);
        continue;
      }
      plane.fill<MuonSpacePoint>(*sp, axisRanges, phiHoughParam_strip,
                                 phiHoughWidth_strip, sp, sp->id().detLayer());
    }
    const MaximumVec_t phiMaxima = peakFinder.findPeaks(plane, axisRanges);
    for (const Maximum_t& max : phiMaxima) {
      MuonHoughMaximum::HitVec hits{max.hitIdentifiers.begin(),
                                    max.hitIdentifiers.end()};
      hits.insert(hits.end(), etaHits.begin(), etaHits.end());
      outMaxima.emplace_back(max.x, max.y, etaMax.tanBeta(),
                             etaMax.interceptY(), std::move(hits));
    }
    MuonId bucketId{etaMax.hits().front()->id()};
    bucketId.setCoordFlags(false, true);
    displayMaxima(ctx, bucketId, phiMaxima, plane, axisRanges);
  }
  return outMaxima;
}

void MuonHoughSeeder::displayMaxima(const AlgorithmContext& ctx,
                                    const MuonId& bucketId,
                                    const MaximumVec_t& maxima,
                                    const HoughPlane_t& plane,
                                    const AxisRange_t& axis) const {
  if (!m_outCanvas) {
    return;
  }
  static std::mutex canvasMutex{};
  // read the hits and circles
  const MuonSegmentContainer& gotTruthSegs = m_inputTruthSegs(ctx);
  std::lock_guard guard{canvasMutex};
  std::vector<std::unique_ptr<TObject>> primitives;

  /// Save the hough accumulator as histogram
  TH2D houghHistoForPlot("houghHist", "HoughPlane;tan(#alpha);z0 [mm]",
                         plane.nBinsX(), axis.xMin, axis.xMax, plane.nBinsY(),
                         axis.yMin, axis.yMax);
  houghHistoForPlot.SetTitle(Form("Station %s, side %s, sector %2d",
                                  to_string(bucketId.msStation()).c_str(),
                                  to_string(bucketId.side()).c_str(),
                                  bucketId.sector()));

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
  MuonSegmentContainer::const_iterator truthItr = gotTruthSegs.begin();
  const int trueCoord = bucketId.measuresEta() ? Acts::eY : Acts::eX;
  while ((truthItr = std::find_if(truthItr, gotTruthSegs.end(),
                                  [bucketId](const MuonSegment& seg) {
                                    return seg.id().sameStation(bucketId);
                                  })) != gotTruthSegs.end()) {
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
  for (const Maximum_t& max : maxima) {
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
  }
  primitives.emplace_back(std::move(legend));
  for (auto& prim : primitives) {
    prim->Draw();
  }

  m_outCanvas->SaveAs("HoughHistograms.pdf");
}
ProcessCode MuonHoughSeeder::initialize() {
  // book the output canvas
  m_outCanvas = std::make_unique<TCanvas>("canvas", "", 800, 800);
  m_outCanvas->SaveAs("HoughHistograms.pdf[");
  m_outCanvas->SetRightMargin(0.12);
  m_outCanvas->SetLeftMargin(0.12);
  gStyle->SetPalette(kGreyScale);
  gStyle->SetOptStat(0);
  return ProcessCode::SUCCESS;
}
ProcessCode MuonHoughSeeder::finalize() {
  m_outCanvas->SaveAs("HoughHistograms.pdf]");
  return ProcessCode::SUCCESS;
}
}  // namespace ActsExamples

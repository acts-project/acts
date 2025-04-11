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
#include "TLatex.h"
#include "TLegend.h"
#include "TMarker.h"
#include "TStyle.h"

namespace {
std::mutex canvasMutex{};
}

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

using PatternSeed = std::pair<double, double>;  // y0, tan theta

ProcessCode MuonHoughSeeder::execute(const AlgorithmContext& ctx) const {
  // read the hits and circles
  const MuonSegmentContainer& gotTruthSegs = m_inputTruthSegs(ctx);
  const MuonSpacePointContainer& gotSpacePoints = m_inputSpacePoints(ctx);

  // configure the binning of the hough plane
  Acts::HoughTransformUtils::HoughPlaneConfig etaPlaneCfg;
  etaPlaneCfg.nBinsX = 25;
  etaPlaneCfg.nBinsY = 25;

  Acts::HoughTransformUtils::HoughPlaneConfig phiPlaneCfg;
  phiPlaneCfg.nBinsX = 10;
  phiPlaneCfg.nBinsY = 10;
  // Acts

  // instantiate the peak finder

  Acts::HoughTransformUtils::PeakFinders::IslandsAroundMaxConfig
      phiPeakFinderCfg;
  phiPeakFinderCfg.fractionCutoff = 0.7;
  phiPeakFinderCfg.threshold = 2.;
  phiPeakFinderCfg.minSpacingBetweenPeaks = {0., 30.};

  // instantiate the hough plane
  HoughPlane_t etaPlane{etaPlaneCfg}, phiPlane{phiPlaneCfg};
  // also instantiate the peak finder
  // PeakFinder_t etaPeakFinder{peakFinderCfg}, phiPeakFinder{phiPeakFinderCfg};

  MuonHoughMaxContainer outMaxima{};

  for (const MuonSpacePointContainer::value_type& bucket : gotSpacePoints) {
    MuonHoughMaxContainer etaMaxima  = constructEtaMaxima(ctx, bucket, etaPlane);
  /*  etaPlane.reset();
    /// First loop over the bucket to find the covered range in x & y
    etaAxisRanges.yMin = 100. * Acts::UnitConstants::m;
    etaAxisRanges.yMax = -100. * Acts::UnitConstants::m;
    phiAxisRanges.yMin = 100. * Acts::UnitConstants::m;
    phiAxisRanges.yMax = -100. * Acts::UnitConstants::m;
    for (const MuonSpacePoint& sp : bucket) {
      if (sp.id().measuresEta()) {
        etaAxisRanges.yMin =
            std::min(etaAxisRanges.yMin,
                     sp.localPosition().y() - m_cfg.etaPlaneMarginIcept);
        etaAxisRanges.yMax =
            std::max(etaAxisRanges.yMax,
                     sp.localPosition().y() + m_cfg.etaPlaneMarginIcept);
      }
      if (sp.id().measuresPhi()) {
        phiAxisRanges.yMin =
            std::min(phiAxisRanges.yMin,
                     sp.localPosition().x() - m_cfg.phiPlaneMarginIcept);
        phiAxisRanges.yMax =
            std::max(phiAxisRanges.yMax,
                     sp.localPosition().x() + m_cfg.phiPlaneMarginIcept);
      }
    }

    for (const MuonSpacePoint& sp : bucket) {
      if (sp.id().technology() == MuonSpacePoint::MuonId::TechField::Mdt) {
        etaPlane.fill<MuonSpacePoint>(sp, etaAxisRanges, etaHoughParamDC_left,
                                      houghWidth_fromDC, &sp,
                                      sp.id().detLayer());
        etaPlane.fill<MuonSpacePoint>(sp, etaAxisRanges, etaHoughParamDC_right,
                                      houghWidth_fromDC, &sp,
                                      sp.id().detLayer());
      } else if (sp.id().measuresEta()) {
        etaPlane.fill<MuonSpacePoint>(sp, etaAxisRanges, etaHoughParam_strip,
                                      etaHoughWidth_strip, &sp,
                                      sp.id().detLayer());
      }
    }
    const auto bucketId{bucket.front().id()};
    m_outCanvas->SetTitle(Form("Station %s, side %s, sector %2d",
                               to_string(bucketId.msStation()).c_str(),
                               to_string(bucketId.side()).c_str(),
                               bucketId.sector()));

    /// Save the hough accumulator as histogram
    TH2D houghHistoForPlot("houghHist", "HoughPlane;tan(#alpha);z0 [mm]",
                           etaPlane.nBinsX(), etaAxisRanges.xMin,
                           etaAxisRanges.xMax, etaPlane.nBinsY(),
                           etaAxisRanges.yMin, etaAxisRanges.yMax);
    houghHistoForPlot.SetTitle(Form("Station %s, side %s, sector %2d",
                                    to_string(bucketId.msStation()).c_str(),
                                    to_string(bucketId.side()).c_str(),
                                    bucketId.sector()));

    for (int bx = 0; bx < houghHistoForPlot.GetNbinsX(); ++bx) {
      for (int by = 0; by < houghHistoForPlot.GetNbinsY(); ++by) {
        houghHistoForPlot.SetBinContent(bx + 1, by + 1, etaPlane.nHits(bx, by));
      }
    }
    /// Set the contours
    int maxHitsAsInt = static_cast<int>(etaPlane.maxHits());
    houghHistoForPlot.SetContour(maxHitsAsInt + 1);
    for (int k = 0; k < maxHitsAsInt + 1; ++k) {
      houghHistoForPlot.SetContourLevel(k, k - 0.5);
    }
    /// Fetch the first set of peaks
    const MaximumVec_t maxima =
        etaPeakFinder.findPeaks(etaPlane, etaAxisRanges);

    std::vector<std::unique_ptr<TObject>> primitives;

    houghHistoForPlot.Draw("COLZ");

    auto legend = std::make_unique<TLegend>(
        0.5, 0.7, 1. - gPad->GetRightMargin(), 1. - gPad->GetTopMargin());
    legend->SetBorderSize(0);
    legend->SetFillStyle(0);
    /// Fetch the truth parameters
    MuonSegmentContainer::const_iterator truthItr = gotTruthSegs.begin();
    while ((truthItr = std::find_if(truthItr, gotTruthSegs.end(),
                                    [bucketId](const MuonSegment& seg) {
                                      return seg.id().sameStation(bucketId);
                                    })) != gotTruthSegs.end()) {
      const MuonSegment& truthSeg{*truthItr};
      const float tanAlpha =
          truthSeg.localDirection().y() / truthSeg.localDirection().z();
      const float intercept = truthSeg.localPosition().y();

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
    {
      auto tl = std::make_unique<TLatex>(
          gPad->GetLeftMargin() + 0.03, 1. - gPad->GetTopMargin() - 0.1,
          Form("Space points in station %lu", bucket.size()));

      tl->SetTextFont(43);
      tl->SetTextSize(24);
      tl->SetNDC();
      primitives.emplace_back(std::move(tl));
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
      primitives.push_back(std::move(marker));
      primitives.push_back(std::move(box));

      MuonHoughMaximum::HitVec hits{max.hitIdentifiers.begin(),
                                    max.hitIdentifiers.end()};
      /// Add all the space points that don't measure the precision coordinate
      for (const MuonSpacePoint& sp : bucket) {
        if (!sp.id().measuresEta()) {
          hits.push_back(&sp);
        }
      }
      etaMaxima.emplace_back(max.x, max.y, std::move(hits));
    }
    primitives.emplace_back(std::move(legend));
    for (auto& prim : primitives) {
      prim->Draw();
    }

    m_outCanvas->SaveAs("HoughHistograms.pdf");
    /// Attempt to find the parameters in the non-precision plane
    for (MuonHoughMaximum& etaMax : etaMaxima) {
      phiPlane.reset();
      unsigned nPhiHits{0};
      for (const MuonSpacePoint* sp : etaMax.hits()) {
        if (sp->id().measuresPhi()) {
          phiPlane.fill<MuonSpacePoint>(*sp, phiAxisRanges, phiHoughParam_strip,
                                        phiHoughWidth_strip, sp,
                                        sp->id().detLayer());
          ++nPhiHits;
        }
      }
      // 2 hits always share a common hough maximum
      if (nPhiHits < 2) {
        outMaxima.push_back(std::move(etaMax));
        continue;
      }
      /// Fetch the first set of peaks
      const MaximumVec_t phiMaxima =
          phiPeakFinder.findPeaks(phiPlane, phiAxisRanges);
      for (const Maximum_t& max : phiMaxima) {
        MuonHoughMaximum::HitVec hits{max.hitIdentifiers.begin(),
                                      max.hitIdentifiers.end()};
        /// Copy all hits not constraining the phi coordinate
        std::ranges::copy_if(
            etaMax.hits(), std::back_inserter(hits),
            [](const MuonSpacePoint* sp) { return !sp->id().measuresPhi(); });
        outMaxima.emplace_back(max.x, max.y, etaMax.tanBeta(),
                               etaMax.interceptY(), std::move(hits));
      }
    } */
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
  /** @brief Addapt the axis ranges of the y intercept */
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
                                    etaHoughWidth_strip, &sp,
                                    sp.id().detLayer());
    }
  }
  /** Extract the maxima from the peak */
  PeakFinderCfg_t peakFinderCfg{};
  peakFinderCfg.fractionCutoff = 0.7;
  peakFinderCfg.threshold = 3.;
  peakFinderCfg.minSpacingBetweenPeaks = {0., 30.};

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
void MuonHoughSeeder::displayMaxima(const AlgorithmContext& ctx,
                                   const MuonId& bucketId,
                                   const MaximumVec_t& maxima,
                                   const HoughPlane_t& plane,
                                   const AxisRange_t& axis) const {
  if (!m_outCanvas) {
    return;
  }
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
  int maxHitsAsInt = static_cast<int>(plane.maxHits());
  houghHistoForPlot.SetContour(maxHitsAsInt + 1);
  for (int k = 0; k < maxHitsAsInt + 1; ++k) {
    houghHistoForPlot.SetContourLevel(k, k - 0.5);
  }

  auto legend = std::make_unique<TLegend>(
    0.5, 0.7, 1. - gPad->GetRightMargin(), 1. - gPad->GetTopMargin());
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
  const float tanAlpha = truthSeg.localDirection()[trueCoord] / truthSeg.localDirection().z();
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
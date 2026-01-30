// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/Io/Root/RootTrackFitterPerformanceWriter.hpp"

#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/Utilities/Helpers.hpp"
#include "ActsExamples/EventData/SimParticle.hpp"
#include "ActsExamples/EventData/Track.hpp"
#include "ActsExamples/Framework/AlgorithmContext.hpp"
#include "ActsExamples/Validation/TrackClassification.hpp"
#include "ActsFatras/EventData/Barcode.hpp"
#include "ActsPlugins/Root/HistogramConverter.hpp"

#include <cstddef>
#include <ostream>
#include <stdexcept>
#include <utility>
#include <vector>

#include <TEfficiency.h>
#include <TFile.h>
#include <TFitResult.h>
#include <TFitResultPtr.h>
#include <TH1.h>
#include <TH2.h>
#include <TProfile.h>

using Acts::VectorHelpers::eta;
using Acts::VectorHelpers::phi;

ActsExamples::RootTrackFitterPerformanceWriter::
    RootTrackFitterPerformanceWriter(
        ActsExamples::RootTrackFitterPerformanceWriter::Config config,
        Acts::Logging::Level level)
    : WriterT(config.inputTracks, "RootTrackFitterPerformanceWriter", level),
      m_cfg(std::move(config)),
      m_resPlotTool(m_cfg.resPlotToolConfig, level),
      m_effPlotTool(m_cfg.effPlotToolConfig, level),
      m_trackSummaryPlotTool(m_cfg.trackSummaryPlotToolConfig, level) {
  // trajectories collection name is already checked by base ctor
  if (m_cfg.inputParticles.empty()) {
    throw std::invalid_argument("Missing particles input collection");
  }
  if (m_cfg.inputTrackParticleMatching.empty()) {
    throw std::invalid_argument("Missing input track particles matching");
  }
  if (m_cfg.filePath.empty()) {
    throw std::invalid_argument("Missing output filename");
  }

  m_inputParticles.initialize(m_cfg.inputParticles);
  m_inputTrackParticleMatching.initialize(m_cfg.inputTrackParticleMatching);

  // the output file can not be given externally since TFile accesses to the
  // same file from multiple threads are unsafe.
  // must always be opened internally
  auto path = m_cfg.filePath;
  m_outputFile = TFile::Open(path.c_str(), "RECREATE");
  if (m_outputFile == nullptr) {
    throw std::invalid_argument("Could not open '" + path + "'");
  }
}

ActsExamples::RootTrackFitterPerformanceWriter::
    ~RootTrackFitterPerformanceWriter() {
  delete m_outputFile;
}

ActsExamples::ProcessCode
ActsExamples::RootTrackFitterPerformanceWriter::finalize() {
  if (m_outputFile == nullptr) {
    return ProcessCode::SUCCESS;
  }

  m_outputFile->cd();

  // Helper lambda to write 2D histogram and extract mean/width histograms
  auto writeWithRefinement = [](TH2F* hist2d, const std::string& meanPrefix,
                                const std::string& widthPrefix) {
    hist2d->Write();

    // Get the histogram name and extract the suffix (e.g., "_d0_vs_eta")
    std::string baseName = hist2d->GetName();
    std::string suffix = baseName.substr(baseName.find('_'));

    // Create mean and width histograms with same X binning as the 2D histogram
    int nBinsX = hist2d->GetNbinsX();
    auto meanHist =
        std::make_unique<TH1F>((meanPrefix + suffix).c_str(),
                 (std::string(hist2d->GetTitle()) + " mean").c_str(), nBinsX,
                 hist2d->GetXaxis()->GetXmin(), hist2d->GetXaxis()->GetXmax());
    auto widthHist =
        std::make_unique<TH1F>((widthPrefix + suffix).c_str(),
                 (std::string(hist2d->GetTitle()) + " width").c_str(), nBinsX,
                 hist2d->GetXaxis()->GetXmin(), hist2d->GetXaxis()->GetXmax());

    // Copy X axis bin edges for variable binning
    if (hist2d->GetXaxis()->GetXbins()->GetSize() > 0) {
      meanHist->SetBins(nBinsX, hist2d->GetXaxis()->GetXbins()->GetArray());
      widthHist->SetBins(nBinsX, hist2d->GetXaxis()->GetXbins()->GetArray());
    }

    // Project each X bin and extract mean/width via Gaussian fit
    for (int j = 1; j <= nBinsX; j++) {
      TH1D* proj = hist2d->ProjectionY(
          Form("%s_projy_bin%d", baseName.c_str(), j), j, j);
      if (proj->GetEntries() > 0) {
        TFitResultPtr r = proj->Fit("gaus", "QS0");
        if ((r.Get() != nullptr) && ((r->Status() % 1000) == 0)) {
          // fill the mean and width into 'j'th bin of the meanHist and
          // widthHist, respectively
          meanHist->SetBinContent(j, r->Parameter(1));
          meanHist->SetBinError(j, r->ParError(1));
          widthHist->SetBinContent(j, r->Parameter(2));
          widthHist->SetBinError(j, r->ParError(2));
        }
      }
      delete proj;
    }

    meanHist->Write();
    widthHist->Write();
    delete meanHist;
    delete widthHist;
  };

  // Write residual histograms
  for (const auto& [name, hist] : m_resPlotTool.res()) {
    ActsPlugins::toRoot(hist)->Write();
  }
  for (const auto& [name, hist] : m_resPlotTool.resVsEta()) {
    writeWithRefinement(ActsPlugins::toRoot(hist), "resmean", "reswidth");
  }
  for (const auto& [name, hist] : m_resPlotTool.resVsPt()) {
    writeWithRefinement(ActsPlugins::toRoot(hist), "resmean", "reswidth");
  }

  // Write pull histograms
  for (const auto& [name, hist] : m_resPlotTool.pull()) {
    ActsPlugins::toRoot(hist)->Write();
  }
  for (const auto& [name, hist] : m_resPlotTool.pullVsEta()) {
    writeWithRefinement(ActsPlugins::toRoot(hist), "pullmean", "pullwidth");
  }
  for (const auto& [name, hist] : m_resPlotTool.pullVsPt()) {
    writeWithRefinement(ActsPlugins::toRoot(hist), "pullmean", "pullwidth");
  }

  // Write efficiency histograms
  ActsPlugins::toRoot(m_effPlotTool.trackEffVsEta())->Write();
  ActsPlugins::toRoot(m_effPlotTool.trackEffVsPhi())->Write();
  ActsPlugins::toRoot(m_effPlotTool.trackEffVsPt())->Write();
  ActsPlugins::toRoot(m_effPlotTool.trackEffVsLogPt())->Write();
  ActsPlugins::toRoot(m_effPlotTool.trackEffVsLowPt())->Write();
  ActsPlugins::toRoot(m_effPlotTool.trackEffVsD0())->Write();
  ActsPlugins::toRoot(m_effPlotTool.trackEffVsZ0())->Write();
  ActsPlugins::toRoot(m_effPlotTool.trackEffVsDeltaR())->Write();
  ActsPlugins::toRoot(m_effPlotTool.trackEffVsProdR())->Write();
  ActsPlugins::toRoot(m_effPlotTool.trackEffVsEtaPhi())->Write();
  ActsPlugins::toRoot(m_effPlotTool.trackEffVsEtaPt())->Write();

  for (const auto& eff : m_effPlotTool.trackEffVsEtaInPtRanges()) {
    ActsPlugins::toRoot(eff)->Write();
  }
  for (const auto& eff : m_effPlotTool.trackEffVsPtInAbsEtaRanges()) {
    ActsPlugins::toRoot(eff)->Write();
  }

  // Write track summary histograms
  ActsPlugins::toRoot(m_trackSummaryPlotTool.nStatesVsEta())->Write();
  ActsPlugins::toRoot(m_trackSummaryPlotTool.nMeasurementsVsEta())->Write();
  ActsPlugins::toRoot(m_trackSummaryPlotTool.nHolesVsEta())->Write();
  ActsPlugins::toRoot(m_trackSummaryPlotTool.nOutliersVsEta())->Write();
  ActsPlugins::toRoot(m_trackSummaryPlotTool.nSharedHitsVsEta())->Write();
  ActsPlugins::toRoot(m_trackSummaryPlotTool.nStatesVsPt())->Write();
  ActsPlugins::toRoot(m_trackSummaryPlotTool.nMeasurementsVsPt())->Write();
  ActsPlugins::toRoot(m_trackSummaryPlotTool.nHolesVsPt())->Write();
  ActsPlugins::toRoot(m_trackSummaryPlotTool.nOutliersVsPt())->Write();
  ActsPlugins::toRoot(m_trackSummaryPlotTool.nSharedHitsVsPt())->Write();

  ACTS_INFO("Wrote performance plots to '" << m_outputFile->GetPath() << "'");

  if (m_outputFile != nullptr) {
    m_outputFile->Close();
  }
  return ProcessCode::SUCCESS;
}

ActsExamples::ProcessCode
ActsExamples::RootTrackFitterPerformanceWriter::writeT(
    const AlgorithmContext& ctx, const ConstTrackContainer& tracks) {
  // Read truth input collections
  const auto& particles = m_inputParticles(ctx);
  const auto& trackParticleMatching = m_inputTrackParticleMatching(ctx);

  // Truth particles with corresponding reconstructed tracks
  std::vector<ActsFatras::Barcode> reconParticleIds;
  reconParticleIds.reserve(particles.size());
  // For each particle within a track, how many hits did it contribute
  std::vector<ParticleHitCount> particleHitCounts;

  // Exclusive access to the tree while writing
  std::lock_guard<std::mutex> lock(m_writeMutex);

  // Loop over all tracks
  for (const auto& track : tracks) {
    // Select reco track with fitted parameters
    if (!track.hasReferenceSurface()) {
      ACTS_WARNING("No fitted track parameters.");
      continue;
    }
    Acts::BoundTrackParameters fittedParameters =
        track.createParametersAtReference();

    // Get the truth-matched particle
    auto imatched = trackParticleMatching.find(track.index());
    if (imatched == trackParticleMatching.end()) {
      ACTS_DEBUG("No truth particle associated with this track, index = "
                 << track.index() << " tip index = " << track.tipIndex());
      continue;
    }
    const auto& particleMatch = imatched->second;

    if (!particleMatch.particle.has_value()) {
      ACTS_DEBUG("No truth particle associated with this track.");
      continue;
    }

    // Get the barcode of the majority truth particle
    SimBarcode majorityParticleId = particleMatch.particle.value();

    // Find the truth particle via the barcode
    auto ip = particles.find(majorityParticleId);
    if (ip == particles.end()) {
      ACTS_DEBUG("Majority particle not found in the particles collection.");
      continue;
    }

    // Record this majority particle ID of this trajectory
    reconParticleIds.push_back(ip->particleId());
    // Fill the residual plots
    m_resPlotTool.fill(ctx.geoContext, ip->initialState(), fittedParameters);
    // Fill the trajectory summary info
    m_trackSummaryPlotTool.fill(fittedParameters, track.nTrackStates(),
                                track.nMeasurements(), track.nOutliers(),
                                track.nHoles(), track.nSharedHits());
  }

  // Fill the efficiency, defined as the ratio between number of tracks with
  // fitted parameter and total truth tracks (assumes one truth particle has
  // one truth track)
  for (const auto& particle : particles) {
    bool isReconstructed = false;
    // Check if the particle has been reconstructed
    if (Acts::rangeContainsValue(reconParticleIds, particle.particleId())) {
      isReconstructed = true;
    }
    // Loop over all the other truth particle and find the distance to the
    // closest one
    double minDeltaR = -1;
    for (const auto& closeParticle : particles) {
      if (closeParticle.particleId() == particle.particleId()) {
        continue;
      }
      double p_phi = phi(particle.direction());
      double p_eta = eta(particle.direction());
      double c_phi = phi(closeParticle.direction());
      double c_eta = eta(closeParticle.direction());
      double distance = sqrt(pow(p_phi - c_phi, 2) + pow(p_eta - c_eta, 2));
      if (minDeltaR == -1 || distance < minDeltaR) {
        minDeltaR = distance;
      }
    }
    m_effPlotTool.fill(ctx.geoContext, particle.initialState(), minDeltaR,
                       isReconstructed);
  }

  return ProcessCode::SUCCESS;
}

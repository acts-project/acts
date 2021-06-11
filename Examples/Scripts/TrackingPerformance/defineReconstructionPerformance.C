// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <array>
#include <iostream>
#include <string>
#include <tuple>
#include <vector>

#include <TCanvas.h>
#include <TEfficiency.h>
#include <TFile.h>
#include <TProfile.h>

#include "TreeAccessors.h"
#include "CommonUtils.h"

/// This script/function reads all the reconstructed tracks from e.g.
/// 'tracksummary_ckf.root' and the (possibly selected) truth particles from
/// e.g. 'track_finder_particles.root' (which contains the info of 'nHits'), and
/// defines the efficiency, fake rate and duplicaiton rate. It aims to make
/// custom definition and tuning of the reconstruction performance easier.
///
/// NB: It's very likely that fiducal cuts are already imposed on the truth
/// particles. Please check the selection criteria in the truth fitting example
/// which writes out the 'track_finder_particles.root'. For instance, if the
/// truth particles are already required to have pT > 1 GeV, it does not make
/// sense to have ptMin = 0.5 GeV here.
///
void defineReconstructionPerformance(
    const std::string& inputTrackSummaryFileName = "tracksummary_ckf.root",
    const std::string& inputSimParticleFileName = "performance_track_finder.root",
    const std::string& trackSummaryTreeName = "tracksummary_ckf",
    const std::string& simParticleTreeName = "track_finder_particles",
    unsigned int nHitsMin = 9, unsigned int nMeasurementsMin = 6,
    double ptMin = 0.5, double truthMatchProbMin = 0.5) {
  // Open the files for reading
  TFile* trackFile = TFile::Open(inputTrackSummaryFileName.c_str(), "read");
  TFile* particleFile = TFile::Open(inputSimParticleFileName.c_str(), "read");

  // Define variables for tree reading (setBranchAddress etc. is done in the
  // contructor)
  TrackVars tVars((TTree*)trackFile->Get(trackSummaryTreeName.c_str()));
  ParticleVars pVars((TTree*)particleFile->Get(simParticleTreeName.c_str()));

  size_t nEvents = tVars.tree->GetEntries();
  std::cout << "There are " << nEvents << " events in total." << std::endl;

  // Define the efficiency plots
  TEfficiency* trackEff_vs_eta = new TEfficiency(
      "trackeff_vs_eta", "Tracking efficiency;Truth #eta [GeV/c];Efficiency", 40, -4, 4);
  TEfficiency* fakeRate_vs_eta = new TEfficiency(
      "fakerate_vs_eta", "Fake rate;#eta [GeV/c];fake rate", 40, -4, 4);
  TEfficiency* duplicateRate_vs_eta = new TEfficiency(
      "duplicaterate_vs_eta", "Duplicate rate;#eta [GeV/c];Duplicate rate", 40, -4, 4);
  TEfficiency* trackEff_vs_pt = new TEfficiency(
      "trackeff_vs_pt", "Tracking efficiency;Truth pt [GeV/c];Efficiency", 40, 0, 100);
  TEfficiency* fakeRate_vs_pt = new TEfficiency(
      "fakerate_vs_pt", "Fake rate;pt [GeV/c];fake rate", 40, 0, 100);
  TEfficiency* duplicateRate_vs_pt = new TEfficiency(
      "duplicaterate_vs_pt", "Duplicate rate;pt [GeV/c];Duplicate rate", 40, 0, 100);

  // Set styles
  setEffStyle(trackEff_vs_eta, 1);
  setEffStyle(fakeRate_vs_eta, 1);
  setEffStyle(duplicateRate_vs_eta, 1);
  setEffStyle(trackEff_vs_pt, 1);
  setEffStyle(fakeRate_vs_pt, 1);
  setEffStyle(duplicateRate_vs_pt, 1);

  // The container from track-particle matching info (Flushed per event)
  std::map<uint64_t, std::vector<RecoTrackInfo>> matchedParticles;

  // Get the sorted entry numbers
  const auto& tEntryNumbers = tVars.entryNumbers;

  // Loop over the events to fill plots
  for (size_t i = 0; i < nEvents; ++i) {
    if (i % 10 == 0) {
      std::cout << "Processed events: " << i << std::endl;
    }

    // Get the tracks
    auto tEntry = tEntryNumbers[i];
    tVars.tree->GetEvent(tEntry);

    // Get the particles
    auto particles = pVars.getParticles(i);
    //    std::cout << "There are " << particles.size()
    //              << " truth particles in event: " << i << std::endl;

    // Loop over the tracks
    // The fake rate is defined as the ratio of selected truth-matched tracks
    // over all selected tracks
    for (size_t j = 0; j < tVars.nStates->size(); ++j) {
      bool hasFittedParameters = tVars.hasFittedParams->at(j);
      auto nMeasurements = tVars.nMeasurements->at(j);
      auto nOutliers = tVars.nOutliers->at(j);
      auto nHoles = tVars.nHoles->at(j);
      auto theta = tVars.eTHETA_fit->at(j);
      auto qop = tVars.eQOP_fit->at(j);
      auto pt = std::abs(1 / qop) * std::sin(theta);
      auto eta = std::atanh(std::cos(theta));
      auto nMajorityHits = tVars.nMajorityHits->at(j);
      auto majorityParticleId = tVars.majorityParticleId->at(j);

      // Select the tVars, e.g. you might also want to add cuts on the
      // nOutliers, nHoles
      if ((!hasFittedParameters) or nMeasurements < nMeasurementsMin or
          pt < ptMin) {
        continue;
      }

      // Fill the fake rate plots
      if (nMajorityHits * 1. / nMeasurements >= truthMatchProbMin) {
        matchedParticles[majorityParticleId].push_back(
            {eta, pt, nMajorityHits, nMeasurements});
        fakeRate_vs_eta->Fill(false, eta);
        fakeRate_vs_pt->Fill(false, pt);
      } else {
        fakeRate_vs_eta->Fill(true, eta);
        fakeRate_vs_pt->Fill(true, pt);
      }
    }  // end of all tracks

    // Loop over all selected and truth-matched tracks
    // The duplicate rate is defined as the ratio of duplicate tracks among all
    // the selected truth-matched tracks (only one track is 'real'; others are
    // 'duplicated')
    for (auto& [id, matchedTracks] : matchedParticles) {
      // Sort all tracks matched to this particle according to majority prob and
      // track quality
      std::sort(matchedTracks.begin(), matchedTracks.end(),
                [](const RecoTrackInfo& lhs, const RecoTrackInfo& rhs) {
                  if (lhs.nMajorityHits > rhs.nMajorityHits) {
                    return true;
                  }
                  if (lhs.nMajorityHits < rhs.nMajorityHits) {
                    return false;
                  }
                  if (lhs.nMeasurements > rhs.nMeasurements) {
                    return true;
                  }
                  return false;
                });
      // Fill the duplication rate plots
      for (size_t k = 0; k < matchedTracks.size(); ++k) {
        auto eta = matchedTracks[k].eta;
        auto pt = matchedTracks[k].pt;
        if (k == 0) {
          duplicateRate_vs_eta->Fill(true, eta);
          duplicateRate_vs_pt->Fill(true, pt);
        } else {
          duplicateRate_vs_eta->Fill(false, eta);
          duplicateRate_vs_pt->Fill(false, pt);
        }
      }
    }  // end of all selected truth-matched tracks

    // Loop over all truth particles in this event
    // The effiency is define as the ratio of selected particles that have been
    // matched with reco
    for (const auto& particle : particles) {
      auto nHits = particle.nHits;
      auto eta = particle.eta;
      auto pt = particle.pt;
      if (nHits < nHitsMin or pt < ptMin) {
        continue;
      }
      uint64_t id = particle.particleId;

      // Fill the efficiency plots
      auto ip = matchedParticles.find(id);
      if (ip != matchedParticles.end()) {
        trackEff_vs_eta->Fill(true, eta);
        trackEff_vs_pt->Fill(true, pt);
      } else {
        trackEff_vs_eta->Fill(false, eta);
        trackEff_vs_pt->Fill(false, pt);
      }
    }  // end of all particles

    matchedParticles.clear();
  }  // end of all events
  std::cout<<"All good. Now plotting..."<< std::endl;

  // Now draw the plots
  TCanvas* c1 = new TCanvas("recoPerf", " ", 1500, 800);
  c1->Divide(3, 2);
  c1->cd(1);
  trackEff_vs_eta->Draw();

  c1->cd(2);
  fakeRate_vs_eta->Draw();

  c1->cd(3);
  duplicateRate_vs_eta->Draw();

  c1->cd(4);
  trackEff_vs_pt->Draw();

  c1->cd(5);
  fakeRate_vs_pt->Draw();

  c1->cd(6);
  duplicateRate_vs_pt->Draw();
}

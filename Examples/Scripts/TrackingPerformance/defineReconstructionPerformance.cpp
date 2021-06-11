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

#include "CommonUtils.h"

struct TreeVars {
  TTree* tree = nullptr;

  /// Event identifier.
  uint32_t eventId;

  /// The entry numbers for accessing events in increased order (there could be
  /// multiple entries corresponding to one event number)
  std::vector<long long> entryNumbers = {};

  virtual void setBranch() = 0;

  virtual void setEntryNumbers() = 0;
};

struct TrackVars : public TreeVars {
  std::vector<unsigned int>* nStates = new std::vector<unsigned int>;
  std::vector<unsigned int>* nMeasurements = new std::vector<unsigned int>;
  std::vector<unsigned int>* nOutliers = new std::vector<unsigned int>;
  std::vector<unsigned int>* nHoles = new std::vector<unsigned int>;
  std::vector<float>* chi2Sum = new std::vector<float>;
  std::vector<std::vector<double>>* chi2OnMeasurements = new std::vector<std::vector<double>>;
  std::vector<unsigned int>* NDF = new std::vector<unsigned int>;
  std::vector<unsigned int>* nMajorityHits = new std::vector<unsigned int>;
  std::vector<uint64_t>* majorityParticleId = new std::vector<uint64_t>;

  std::vector<bool>* hasFittedParams = new std::vector<bool>;
  std::vector<float>* eLOC0_fit = new std::vector<float>;
  std::vector<float>* eLOC1_fit = new std::vector<float>;
  std::vector<float>* ePHI_fit = new std::vector<float>;
  std::vector<float>* eTHETA_fit = new std::vector<float>;
  std::vector<float>* eQOP_fit = new std::vector<float>;
  std::vector<float>* eT_fit = new std::vector<float>;
  std::vector<float>* err_eLOC0_fit = new std::vector<float>;
  std::vector<float>* err_eLOC1_fit = new std::vector<float>;
  std::vector<float>* err_ePHI_fit = new std::vector<float>;
  std::vector<float>* err_eTHETA_fit = new std::vector<float>;
  std::vector<float>* err_eQOP_fit = new std::vector<float>;
  std::vector<float>* err_eT_fit = new std::vector<float>;

  void setBranch() {
    tree->SetBranchAddress("event_nr", &eventId);
    tree->SetBranchAddress("nStates", &nStates);
    tree->SetBranchAddress("nMeasurements", &nMeasurements);
    tree->SetBranchAddress("nOutliers", &nOutliers);
    tree->SetBranchAddress("nHoles", &nHoles);
    tree->SetBranchAddress("chi2Sum", &chi2Sum);
    tree->SetBranchAddress("chi2OnMeasurements", &chi2OnMeasurements);
    tree->SetBranchAddress("nMajorityHits", &nMajorityHits);
    tree->SetBranchAddress("majorityParticleId", &majorityParticleId);
    tree->SetBranchAddress("NDF", &NDF);

    tree->SetBranchAddress("hasFittedParams", &hasFittedParams);
    tree->SetBranchAddress("eLOC0_fit", &eLOC0_fit);
    tree->SetBranchAddress("eLOC1_fit", &eLOC1_fit);
    tree->SetBranchAddress("ePHI_fit", &ePHI_fit);
    tree->SetBranchAddress("eTHETA_fit", &eTHETA_fit);
    tree->SetBranchAddress("eQOP_fit", &eQOP_fit);
    tree->SetBranchAddress("eT_fit", &eT_fit);
    tree->SetBranchAddress("err_eLOC0_fit", &err_eLOC0_fit);
    tree->SetBranchAddress("err_eLOC1_fit", &err_eLOC1_fit);
    tree->SetBranchAddress("err_ePHI_fit", &err_ePHI_fit);
    tree->SetBranchAddress("err_eTHETA_fit", &err_eTHETA_fit);
    tree->SetBranchAddress("err_eQOP_fit", &err_eQOP_fit);
    tree->SetBranchAddress("err_eT_fit", &err_eT_fit);
  }

  void setEntryNumbers() {
    entryNumbers.resize(tree->GetEntries());
    tree->Draw("event_nr", "", "goff");
    // Sort to get the entry numbers of the ordered events
    TMath::Sort(tree->GetEntries(), tree->GetV1(), entryNumbers.data(), false);
  }
};

struct ParticleVars : public TreeVars {
  std::vector<uint64_t>* particleId = new std::vector<uint64_t>;
  std::vector<int32_t>* particleType = new std::vector<int32_t>;
  std::vector<uint32_t>* process = new std::vector<uint32_t>;
  std::vector<float>* vx = new std::vector<float>;
  std::vector<float>* vy = new std::vector<float>;
  std::vector<float>* vz = new std::vector<float>;
  std::vector<float>* vt = new std::vector<float>;
  std::vector<float>* px = new std::vector<float>;
  std::vector<float>* py = new std::vector<float>;
  std::vector<float>* pz = new std::vector<float>;
  std::vector<float>* m = new std::vector<float>;
  std::vector<float>* q = new std::vector<float>;
  std::vector<float>* eta = new std::vector<float>;
  std::vector<float>* phi = new std::vector<float>;
  std::vector<float>* pt = new std::vector<float>;
  std::vector<float>* p = new std::vector<float>;
  std::vector<uint32_t>* vertexPrimary = new std::vector<uint32_t>;
  std::vector<uint32_t>* vertexSecondary = new std::vector<uint32_t>;
  std::vector<uint32_t>* particle = new std::vector<uint32_t>;
  std::vector<uint32_t>* generation = new std::vector<uint32_t>;
  std::vector<uint32_t>* subParticle = new std::vector<uint32_t>;
  std::vector<uint32_t>* nHits = new std::vector<uint32_t>;

  void setBranch() {
    // Set the branches
    tree->SetBranchAddress("event_id", &eventId);
    tree->SetBranchAddress("particle_id", &particleId);
    tree->SetBranchAddress("particle_type", &particleType);
    tree->SetBranchAddress("process", &process);
    tree->SetBranchAddress("vx", &vx);
    tree->SetBranchAddress("vy", &vy);
    tree->SetBranchAddress("vz", &vz);
    tree->SetBranchAddress("vt", &vt);
    tree->SetBranchAddress("p", &p);
    tree->SetBranchAddress("px", &px);
    tree->SetBranchAddress("py", &py);
    tree->SetBranchAddress("pz", &pz);
    tree->SetBranchAddress("m", &m);
    tree->SetBranchAddress("q", &q);
    tree->SetBranchAddress("eta", &eta);
    tree->SetBranchAddress("phi", &phi);
    tree->SetBranchAddress("pt", &pt);
    tree->SetBranchAddress("vertex_primary", &vertexPrimary);
    tree->SetBranchAddress("vertex_secondary", &vertexSecondary);
    tree->SetBranchAddress("particle", &particle);
    tree->SetBranchAddress("generation", &generation);
    tree->SetBranchAddress("sub_particle", &subParticle);
    tree->SetBranchAddress("nHits", &nHits);
  }

  void setEntryNumbers() {
    entryNumbers.resize(tree->GetEntries());
    tree->Draw("event_id", "", "goff");
    // Sort to get the entry numbers of the ordered events
    TMath::Sort(tree->GetEntries(), tree->GetV1(), entryNumbers.data(), false);
  }
};

/// This script reads all the reconstructed tracks from e.g.
/// 'tracksummary_ckf.root' and the truth particles from e.g.
/// 'particles_initial.root' and defines the efficiency, fake rate and
/// duplicaiton rate. It aims to make custom definition and tuning of the
/// reconstruction performance easier.
///
void defineReconstructionPerformance(
    const std::string& inputTrackSummaryFileName =
        "tracksummary_ckf.root",
    const std::string& inputTruthParticleFileName =
        "particles_initial.root",
    const std::string& trackSummaryTreeName = "tracksummary_ckf",
    const std::string& truthParticleTreeName = "particles",
    unsigned int nHitsMin = 9, unsigned int nMeasurementsMin = 6,
    double ptMin = 1., double truthMatchProbMin = 0.5) {
  // Open the files for reading
  TFile* trackFile = TFile::Open(inputTrackSummaryFileName.c_str(), "read");
  TFile* particleFile = TFile::Open(inputTruthParticleFileName.c_str(), "read");

  // Define variables for tree reading
  TrackVars tracks;
  ParticleVars particles;
  tracks.tree = (TTree*)trackFile->Get(trackSummaryTreeName.c_str());
  particles.tree = (TTree*)particleFile->Get(truthParticleTreeName.c_str());

  // Check if there are exactly the same number of events in the two files
  size_t nEvents = tracks.tree->GetEntries();
  if (nEvents != particles.tree->GetEntries()) {
    throw std::invalid_argument(
        "The number of events are not consistent in the track summary and "
        "truth particles files");
  }
  std::cout << "There are " << nEvents << " events in total." << std::endl;;

  // Set the branch address for variables
  tracks.setBranch();
  particles.setBranch();

  // Get the entry numbers for events in increasing order (this is to
  // 'synchronize' the events in the two root files)
  tracks.setEntryNumbers();
  particles.setEntryNumbers();

  const auto& tEntryNumbers = tracks.entryNumbers;
  const auto& pEntryNumbers = particles.entryNumbers;

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

  struct RecoTrackInfo {
    double eta;
    double pt;
    unsigned int nMajorityHits;
    unsigned int nMeasurements;
  };

  std::map<uint64_t, std::vector<RecoTrackInfo>> matchedParticles;
  // Loop over the events to fill plots
  for (size_t i = 0; i < nEvents; ++i) {
    if (i % 10 == 0) {
      std::cout << "Processing event: " << i << std::endl;
    }
    auto tEntry = tEntryNumbers[i];
    auto pEntry = pEntryNumbers[i];

    tracks.tree->GetEvent(tEntry);
    particles.tree->GetEvent(pEntry);

    // Loop over the tracks
    // The fake rate is defined as the ratio of selected truth-matched tracks
    // over all selected tracks
    for (size_t j = 0; j < tracks.nStates->size(); ++j) {
      auto nMeasurements = tracks.nMeasurements->at(j);
      auto nOutliers = tracks.nOutliers->at(j);
      auto nHoles = tracks.nHoles->at(j);
      auto theta = tracks.eTHETA_fit->at(j);
      auto qop = tracks.eQOP_fit->at(j);
      auto pt = std::abs(1 / qop) * std::sin(theta);
      auto eta = std::atanh(std::cos(theta));
      auto nMajorityHits = tracks.nMajorityHits->at(j);
      auto majorityParticleId = tracks.majorityParticleId->at(j);

      // Select the tracks, e.g. you might also want to add cuts on the
      // nOutliers, nHoles
      if ( (!tracks.hasFittedParams->at(j)) or nMeasurements < nMeasurementsMin or
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
    for (size_t j = 0; j < particles.particleId->size(); ++j) {
      auto nHits = particles.nHits->at(j);
      auto teta = particles.eta->at(j);
      auto tpt = particles.pt->at(j);
      if (nHits < nHitsMin or tpt < ptMin) {
        continue;
      }

      // Fill the efficiency plots
      auto ip = matchedParticles.find(particles.particleId->at(j));
      if (ip != matchedParticles.end()) {
        trackEff_vs_eta->Fill(true, teta);
        trackEff_vs_pt->Fill(true, tpt);
      } else {
        trackEff_vs_eta->Fill(false, teta);
        trackEff_vs_pt->Fill(false, tpt);
      }
    }  // end of all particles

    matchedParticles.clear();
  }  // end of all events

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

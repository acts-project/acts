// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <array>
#include <string>
#include <vector>

#include <TCanvas.h>
#include <TEfficiency.h>
#include <TFile.h>
#include <TTree.h>


struct RecoTrackInfo {
  double eta;
  double pt;
  unsigned int nMajorityHits;
  unsigned int nMeasurements;
};

struct ParticleInfo {
  ULong64_t particleId = 0;
  double eta = 0;
  double p = 0;
  double pt = 0;
  UShort_t nHits = 0;
};


struct TreeReader {
  // The constructor
  TreeReader(TTree* tree_) : tree(tree_){};

  // The tree being read
  TTree* tree = nullptr;

  /// The entry numbers for accessing events in increased order (there could be
  /// multiple entries corresponding to one event number)
  std::vector<long long> entryNumbers = {};
};


/// Struct used for reading tracks written out by the RootTrajectorySummaryWriter
///
struct TrackReader : public TreeReader {
  // Delete the default constructor
  TrackReader() = delete;

  // The constructor
  TrackReader(TTree* tree_) : TreeReader(tree_) {
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

    // It might not be necessary, but we still do it to be sure
    entryNumbers.resize(tree->GetEntries());
    tree->Draw("event_nr", "", "goff");
    // Sort to get the entry numbers of the ordered events
    TMath::Sort(tree->GetEntries(), tree->GetV1(), entryNumbers.data(), false);
  }

  // The variables
  uint32_t eventId;
  std::vector<unsigned int>* nStates = new std::vector<unsigned int>;
  std::vector<unsigned int>* nMeasurements = new std::vector<unsigned int>;
  std::vector<unsigned int>* nOutliers = new std::vector<unsigned int>;
  std::vector<unsigned int>* nHoles = new std::vector<unsigned int>;
  std::vector<float>* chi2Sum = new std::vector<float>;
  std::vector<std::vector<double>>* chi2OnMeasurements =
      new std::vector<std::vector<double>>;
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
};


/// Struct used for reading particles written out by the TrackFinderPerformanceWriter 
///
struct ParticleReader : public TreeReader {
  // Delete the default constructor
  ParticleReader() = delete;

  // The constructor
  ParticleReader(TTree* tree_) : TreeReader(tree_) {
    tree->SetBranchAddress("event_id", &eventId);
    tree->SetBranchAddress("particle_id", &particleId);
    tree->SetBranchAddress("particle_type", &particleType);
    tree->SetBranchAddress("vx", &vx);
    tree->SetBranchAddress("vy", &vy);
    tree->SetBranchAddress("vz", &vz);
    tree->SetBranchAddress("vt", &vt);
    tree->SetBranchAddress("px", &px);
    tree->SetBranchAddress("py", &py);
    tree->SetBranchAddress("pz", &pz);
    tree->SetBranchAddress("m", &m);
    tree->SetBranchAddress("q", &q);
    tree->SetBranchAddress("nhits", &nHits);
    tree->SetBranchAddress("ntracks", &nTracks);
    tree->SetBranchAddress("ntracks_majority", &nTracksMajority);

    // It might not be necessary, but we still do it to be sure
    entryNumbers.resize(tree->GetEntries());
    tree->Draw("event_id", "", "goff");
    // Sort to get the entry numbers of the ordered events
    TMath::Sort(tree->GetEntries(), tree->GetV1(), entryNumbers.data(), false);
  }

  // Get all the particles with requested event id
  std::vector<ParticleInfo> getParticles(const uint32_t& eventNumber) const {
    // Find the start entry and the batch size for this event
    std::string eventNumberStr = std::to_string(eventNumber);
    std::string findStartEntry = "event_id<" + eventNumberStr;
    std::string findParticlesSize = "event_id==" + eventNumberStr;
    size_t startEntry = tree->GetEntries(findStartEntry.c_str());
    size_t nParticles = tree->GetEntries(findParticlesSize.c_str());
    if (nParticles == 0) {
      throw std::invalid_argument(
          "No particles found. Please check the input file.");
    }
    std::vector<ParticleInfo> particles;
    particles.reserve(nParticles);
    for (unsigned int i = 0; i < nParticles; ++i) {
      tree->GetEvent(entryNumbers[startEntry + i]);
      auto pt = std::hypot(px, py);
      auto p = std::hypot(pt, pz);
      auto eta = std::atanh(pz / p * 1.);
      particles.push_back({particleId, eta, p, pt, nHits});
    }
    return particles;
  }

  // The variables
  ULong64_t eventId = 0;
  ULong64_t particleId = 0;
  Int_t particleType = 0;
  float vx = 0, vy = 0, vz = 0;
  float vt = 0;
  float px = 0, py = 0, pz = 0;
  float m = 0;
  float q = 0;
  UShort_t nHits = 0;
  UShort_t nTracks = 0;
  UShort_t nTracksMajority = 0;
};

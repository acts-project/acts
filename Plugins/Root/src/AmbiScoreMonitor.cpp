// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Plugins/Root/AmbiScoreMonitor.hpp"

#include "Acts/AmbiguityResolution/ScoreBasedAmbiguityResolution.hpp"

#include <TFile.h>
#include <TTree.h>

void Acts::saveScoreMonitor(
    const std::vector<Acts::ScoreBasedAmbiguityResolution::ScoreMonitor>&
        scoreMonitor,
    const std::string& monitorFilePath,
    const std::vector<std::string>& detectorNames,
    const std::vector<std::string>& optionalFunctions) {
  // Open ROOT file for writing
  TFile* file = TFile::Open(monitorFilePath.c_str(), "UPDATE");
  if (!file || file->IsZombie()) {
    throw std::runtime_error("Could not create ROOT file: " + monitorFilePath);
  }

  // Create a TTree
  TTree* tree = new TTree("ScoreMonitorTree", "Score Monitor Data");

  // Variables for branches
  double pT = 0.0;
  double eta = 0.0;
  double phi = 0.0;

  double ptScore = 0;
  std::vector<double> detectorHitScore;
  std::vector<double> detectorHoleScore;
  std::vector<double> detectorOutlierScore;
  std::vector<double> detectorOtherScore;
  double chi2Score = 0;

  std::vector<double> optionalScore;
  std::vector<std::string> detectorNamesroot(detectorNames);
  std::vector<std::string> optionalFunctionsroot(optionalFunctions);

  double totalScore = 0;

  // Create branches
  tree->Branch("pT", &pT);
  tree->Branch("eta", &eta);
  tree->Branch("phi", &phi);
  tree->Branch("ptScore", &ptScore);
  tree->Branch("chi2Score", &chi2Score);
  tree->Branch("totalScore", &totalScore);

  tree->Branch("detectorHitScore", &detectorHitScore);
  tree->Branch("detectorHoleScore", &detectorHoleScore);
  tree->Branch("detectorOutlierScore", &detectorOutlierScore);
  tree->Branch("detectorOtherScore", &detectorOtherScore);
  tree->Branch("optionalScore", &optionalScore);
  tree->Branch("detectorNames", &detectorNamesroot);
  tree->Branch("optionalFunctions", &optionalFunctionsroot);

  // Fill the tree
  for (const auto& monitor : scoreMonitor) {
    pT = monitor.pT;
    eta = monitor.eta;
    phi = monitor.phi;

    ptScore = monitor.ptScore;
    chi2Score = monitor.chi2Score;
    totalScore = monitor.totalScore;

    detectorHitScore = monitor.detectorHitScore;
    detectorHoleScore = monitor.detectorHoleScore;
    detectorOutlierScore = monitor.detectorOutlierScore;
    detectorOtherScore = monitor.detectorOtherScore;
    optionalScore = monitor.optionalScore;
    detectorNamesroot = detectorNames;
    optionalFunctionsroot = optionalFunctions;

    tree->Fill();
  }

  // Write tree to file
  file->Write();
  file->Close();
  delete file;
}

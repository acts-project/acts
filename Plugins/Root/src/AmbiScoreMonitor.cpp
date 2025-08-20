
#include "Acts/Plugins/Root/AmbiScoreMonitor.hpp"

#include "Acts/AmbiguityResolution/ScoreBasedAmbiguityResolution.hpp"
#include <TFile.h>
#include <TTree.h>

void Acts::saveScoreMonitor(
    const std::vector<Acts::ScoreBasedAmbiguityResolution::ScoreMonitor>& scoreMonitor,
    const std::string& monitorFilePath) {
  // Open ROOT file for writing
  TFile* file = TFile::Open(monitorFilePath.c_str(), "RECREATE");
  if (!file || file->IsZombie()) {
    throw std::runtime_error("Could not create ROOT file: " + monitorFilePath);
  }

  // Create a TTree
  TTree* tree = new TTree("ScoreMonitorTree", "Score Monitor Data");

  // Variables for branches
  std::size_t ptScore;
  std::size_t chi2Score;
  std::size_t totalScore;

  std::vector<std::size_t> detectorHitScore;
  std::vector<std::size_t> detectorHoleScore;
  std::vector<std::size_t> detectorOutlierScore;
  std::vector<std::size_t> detectorOtherScore;
  std::vector<std::size_t> optionalScore;

  // Create branches
  tree->Branch("ptScore", &ptScore);
  tree->Branch("chi2Score", &chi2Score);
  tree->Branch("totalScore", &totalScore);

  tree->Branch("detectorHitScore", &detectorHitScore);
  tree->Branch("detectorHoleScore", &detectorHoleScore);
  tree->Branch("detectorOutlierScore", &detectorOutlierScore);
  tree->Branch("detectorOtherScore", &detectorOtherScore);
  tree->Branch("optionalScore", &optionalScore);

  // Fill the tree
  for (const auto& monitor : scoreMonitor) {
    ptScore = monitor.ptScore;
    chi2Score = monitor.chi2Score;
    totalScore = monitor.totalScore;

    detectorHitScore = monitor.detectorHitScore;
    detectorHoleScore = monitor.detectorHoleScore;
    detectorOutlierScore = monitor.detectorOutlierScore;
    detectorOtherScore = monitor.detectorOtherScore;
    optionalScore = monitor.optionalScore;

    tree->Fill();
  }

  // Write tree to file
  file->Write();
  file->Close();
  delete file;
}

// This file is part of the Acts project.
//
// Copyright (C) 2019-2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <array>
#include <iostream>
#include <map>
#include <string>
#include <vector>

#include <TCanvas.h>
#include <TColor.h>
#include <TDirectory.h>
#include <TError.h>
#include <TF1.h>
#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TLegend.h>
#include <TMath.h>
#include <TProfile2D.h>
#include <TStyle.h>
#include <TTree.h>

#include "../CommonUtils.h"

using namespace ROOT;

int vertexPerformance(const std::string& inFile, const std::string& treeName,
                      const std::string& outFile,
                      const std::string& saveAs = "") {
  // Some style options
  gStyle->SetOptFit(0000);
  gStyle->SetOptStat(0000);
  gStyle->SetPadLeftMargin(0.20);
  gStyle->SetPadRightMargin(0.1);
  gStyle->SetPadTopMargin(0.1);
  gStyle->SetPadBottomMargin(0.15);

  // Section 0: file handling ---------------------------------------------
  //
  // Open root file written by RootTrajectoryWriter
  // Create output root file
  std::cout << "Opening file: " << inFile << std::endl;
  TFile* file = TFile::Open(inFile.c_str(), "read");

  // Bail out if no tree was found
  if (file == nullptr) {
    return -1;
  }

  std::cout << "Reading tree: " << treeName << std::endl;
  TTree* tree = (TTree*)file->Get(treeName.c_str());

  // Bail out if no tree was found
  if (tree == nullptr) {
    return -2;
  }

  TFile* out = TFile::Open(outFile.c_str(), "recreate");
  out->cd();

  std::vector<float>* vdiffx{nullptr};
  std::vector<float>* vdiffy{nullptr};
  std::vector<float>* vdiffz{nullptr};

  std::vector<float>* vcovXX{nullptr};

  tree->SetBranchAddress("diffx", &vdiffx);
  tree->SetBranchAddress("diffy", &vdiffy);
  tree->SetBranchAddress("diffz", &vdiffz);
  tree->SetBranchAddress("covXX", &vcovXX);

  auto* h_diffx = new TH1F{"diffx", "diffx", 100, -1, 1};
  auto* h_diffy = new TH1F{"diffy", "diffy", 100, -1, 1};
  auto* h_diffz = new TH1F{"diffz", "diffz", 100, -50, 50};
  auto* h_covXX = new TH1F{"covXX", "covXX", 100, -50, 50};

  auto fill = [](const auto& v, auto* h) {
    for (auto value : v) {
      h->Fill(value);
    }
  };

  for (size_t i = 0; i <= tree->GetEntries(); i++) {
    tree->GetEntry(i);
    // for (auto diffx : *vdiffx) {
    // h_diffx->Fill(diffx);
    // }
    // for (auto diffy : *vdiffy) {
    // h_diffy->Fill(diffy);
    // }
    // for (auto diffy : *vdiffy) {
    // h_diffy->Fill(diffy);
    // }
    fill(*vdiffx, h_diffx);
    fill(*vdiffy, h_diffy);
    fill(*vdiffz, h_diffz);
    fill(*vcovXX, h_covXX);
  }

  h_diffx->Write();
  h_diffy->Write();
  h_diffz->Write();
  h_covXX->Write();

  if (out != nullptr) {
    out->Close();
  }
  return 0;
}

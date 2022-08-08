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
#include <TEfficiency.h>
#include <TTreeReader.h>


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
  /* TTree* tree = (TTree*)file->Get(treeName.c_str()); */

  /* // Bail out if no tree was found */
  /* if (tree == nullptr) { */
    /* return -2; */
  /* } */

  TFile* out = TFile::Open(outFile.c_str(), "recreate");
  out->cd();

  auto* recoOverTrue = new TEfficiency("nRecoOverTrue", "nRecoOverTrue", 1, -0.5, 0.5);

  TTreeReader reader{treeName.c_str(), file};

  TTreeReaderValue<int> nRecoVtx{reader, "nRecoVtx"};

  while(reader.Next()) {
    std::cout << *nRecoVtx << std::endl;
  }


  recoOverTrue->Write();


  if (out != nullptr) {
    out->Close();
  }
  return 0;
}

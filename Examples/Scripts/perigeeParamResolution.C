// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <TF1.h>
#include <TH1F.h>
#include <TMath.h>
#include <TTree.h>
#include <iostream>
#include <map>
#include <vector>

using namespace ROOT;

void
setHistStyle(TH1F* hist, short color);

// This ROOT script will plot the residual and pull of perigee track parameters (d0,
// z0, phi, theta, q/p, t) from root file produced by the TrackFitterPerformanceWriter 
//
void
perigeeParamResolution(const std::string& inFile)
{
  gStyle->SetOptFit(0000);
  gStyle->SetOptStat(0000);
  gStyle->SetPadLeftMargin(0.20);
  gStyle->SetPadRightMargin(0.05);
  gStyle->SetPadTopMargin(0.1);
  gStyle->SetPadBottomMargin(0.15);
 
  // Open root file written by RootTrajectoryWriter
  std::cout << "Opening file: " << inFile << std::endl;
  TFile*  file = TFile::Open(inFile.c_str(), "read");
 
  // Track parameter name
  std::vector<std::string> paramNames
      = {"d0", "z0", "phi", "theta", "qop", "t"};

  map<string, TH1F*> res;
  map<string, TH1F*> pull;

  // Create the hists and set up
  for (const auto& par: paramNames) {
    // residual hists
    res[par] = (TH1F*) file->Get(Form("res_%s", par.c_str()));
    // pull hists
    pull[par] = (TH1F*) file->Get(Form("pull_%s", par.c_str())); 

    // set style
    setHistStyle(res[par], 2);
    setHistStyle(pull[par], 2);
  }


  // plotting residual
  TCanvas* c1 = new TCanvas("c1", "c1", 1200, 800);
  c1->Divide(3, 2);
  for (size_t ipar = 0; ipar < paramNames.size(); ipar++) {
    c1->cd(ipar + 1);
    res[paramNames.at(ipar)]->Draw("");
  }

  // plotting pull
  TCanvas* c2 = new TCanvas("c2", "c2", 1200, 800);
  c2->Divide(3, 2);
  for (size_t ipar = 0; ipar < paramNames.size(); ipar++) {
    c2->cd(ipar + 1);
    pull[paramNames.at(ipar)]->Draw("");
  }
}

// function to set up the histgram style
void
setHistStyle(TH1F* hist, short color = 1)
{
  hist->GetXaxis()->SetTitleSize(0.05);
  hist->GetYaxis()->SetTitleSize(0.05);
  hist->GetXaxis()->SetLabelSize(0.05);
  hist->GetYaxis()->SetLabelSize(0.05);
  hist->GetXaxis()->SetTitleOffset(1.);
  hist->GetYaxis()->SetTitleOffset(1.8);
  hist->GetXaxis()->SetNdivisions(505);
  hist->SetMarkerStyle(20);
  hist->SetMarkerSize(0.8);
  hist->SetLineWidth(2);
  //hist->SetTitle("");
  hist->SetLineColor(1);
  hist->SetMarkerColor(color);
}

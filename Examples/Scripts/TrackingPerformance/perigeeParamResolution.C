// This file is part of the Acts project.
//
// Copyright (C) 2019-2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <iostream>
#include <map>
#include <vector>

#include <TCanvas.h>
#include <TF1.h>
#include <TFile.h>
#include <TH1F.h>
#include <TLegend.h>
#include <TMath.h>
#include <TStyle.h>
#include <TTree.h>

using namespace ROOT;

void setHistStyle(TH1F* hist, short color);

// This ROOT script will plot the residual and pull of perigee track parameters
// (d0, z0, phi, theta, q/p, t) from root file produced by the
// TrackFitterPerformanceWriter
//
int perigeeParamResolution(const std::string& inFile, bool fit = true) {
  gStyle->SetOptFit(0000);
  gStyle->SetOptStat(0000);
  gStyle->SetPadLeftMargin(0.20);
  gStyle->SetPadRightMargin(0.05);
  gStyle->SetPadTopMargin(0.1);
  gStyle->SetPadBottomMargin(0.15);

  // Open root file written by RootTrajectoryWriter
  std::cout << "Opening file: " << inFile << std::endl;
  TFile* file = TFile::Open(inFile.c_str(), "read");

  if (file == nullptr) {
    return -1;
  }

  // Track parameter name
  std::vector<std::array<std::string, 4>> paramNames = {
      {"d0", "(d_{0}^{rec} - d_{0}^{true})", "[mm]", "#sigma(d_{0})"},
      {"z0", "(z_{0}^{rec} - z_{0}^{true})", "[mm]", "#sigma(z_{0})"},
      {"phi", "(#phi^{rec} - #phi^{true})", "", "#sigma(#phi)"},
      {"theta", "(#theta^{rec} - #theta^{true})", "", "#sigma(#theta)"},
      {"qop", "((q/p)^{rec} - (q/p)^{true})", "[GeV^{-1}]", "#sigma(q/p)"},
      {"t", "(t^{rec} - t^{true})", "[ns]", "#sigma(t)"}};

  std::map<std::string, TH1F*> res;
  std::map<std::string, TH1F*> pull;

  // Create the hists and set up
  for (const auto& parterms : paramNames) {
    std::string par = parterms[0];

    // residual hists
    res[par] = (TH1F*)file->Get(Form("res_%s", par.c_str()));
    // pull hists
    pull[par] = (TH1F*)file->Get(Form("pull_%s", par.c_str()));

    // set style
    setHistStyle(res[par], kBlack);
    setHistStyle(pull[par], kBlack);
  }

  // plotting residual
  TCanvas* c1 = new TCanvas("c1", "c1", 1200, 800);
  c1->Divide(3, 2);
  for (size_t ipar = 0; ipar < paramNames.size(); ipar++) {
    c1->cd(ipar + 1);
    auto resHist = res[paramNames.at(ipar)[0]];

    std::string xtitle = paramNames.at(ipar)[1] + std::string(" ");
    xtitle += paramNames.at(ipar)[2];

    Double_t scale = 1. / resHist->Integral();
    resHist->Scale(scale);

    resHist->GetXaxis()->SetTitle(xtitle.c_str());
    resHist->GetYaxis()->SetTitle("Arb. Units");
    resHist->Draw("hist");
  }

  // plotting pull
  TCanvas* c2 = new TCanvas("c2", "c2", 1200, 800);
  c2->Divide(3, 2);
  for (size_t ipar = 0; ipar < paramNames.size(); ipar++) {
    c2->cd(ipar + 1);

    TLegend* legend = new TLegend(0.7, 0.7, 0.9, 0.9);

    auto pullHist = pull[paramNames.at(ipar)[0]];

    Double_t scale = 1. / pullHist->Integral();
    pullHist->Scale(scale);

    if (fit) {
      pullHist->Fit("gaus", "q");
      TF1* gauss = pullHist->GetFunction("gaus");
      float mu = gauss->GetParameter(1);
      float sigma = gauss->GetParameter(2);

      // Draw the sigma
      TString mu_info;
      mu_info += mu;
      TString sigma_info;
      sigma_info += sigma;

      TString mfit_info = "#mu = ";
      mfit_info += mu_info(0, 5);

      TString sfit_info = "#sigma = ";
      sfit_info += sigma_info(0, 5);

      legend->AddEntry(gauss, sfit_info.Data(), "l");
      legend->AddEntry(gauss, mfit_info.Data(), "");
    }

    std::string xtitle = paramNames.at(ipar)[1] + std::string("/");
    xtitle += paramNames.at(ipar)[3];

    pullHist->GetXaxis()->SetTitle(xtitle.c_str());
    pullHist->GetYaxis()->SetTitle("Arb. Units");

    pullHist->Draw("pe");
    legend->SetBorderSize(0);
    legend->SetFillStyle(0);
    legend->SetTextFont(42);
    legend->Draw();
  }
  return 1;
}

// function to set up the histgram style
void setHistStyle(TH1F* hist, short color = 1) {
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
  hist->SetTitle("");
  hist->SetLineColor(1);
  hist->SetMarkerColor(color);
}

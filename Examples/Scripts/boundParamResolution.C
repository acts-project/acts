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

// This ROOT script will plot the residual and pull of track parameters (loc1,
// loc2, phi, theta, q/p, t) from root file produced by the RootTrajectoryWriter
//
void
boundParamResolution(const std::string& inFile,
          const std::string& treeName)
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
  std::cout << "Reading tree: " << treeName << std::endl;
  TTree* tree = (TTree*)file->Get(treeName.c_str());
 
  // Track parameter name
  std::vector<std::string> paramNames
      = {"loc1", "loc2", "#phi", "#theta", "q/p" ,"t"};

  // Residual range
  std::map<std::string, double> paramResidualRange = {{"loc1", 0.5},
                                                      {"loc2", 0.5},
                                                      {"#phi", 0.005},
                                                      {"#theta", 0.005},
                                                      {"q/p", 0.05},
                                                      {"t", 0.005}};
  // Pull range
  double pullRange = 7;

  map<string, TH1F*> res_prt;
  map<string, TH1F*> res_flt;
  map<string, TH1F*> res_smt;
  map<string, TH1F*> pull_prt;
  map<string, TH1F*> pull_flt;
  map<string, TH1F*> pull_smt;

  // Create the hists and set up
  for (const auto& [par, resRange] : paramResidualRange) {
    // residual hists
    res_prt[par] = new TH1F(Form("res_prt_%s", par.c_str()),
                            Form("residual of %s", par.c_str()),
                            100,
                            -1 * resRange,
                            resRange);
    res_flt[par] = new TH1F(Form("res_flt_%s", par.c_str()),
                            Form("residual of %s", par.c_str()),
                            100,
                            -1 * resRange,
                            resRange);
    res_smt[par] = new TH1F(Form("res_smt_%s", par.c_str()),
                            Form("residual of %s", par.c_str()),
                            100,
                            -1 * resRange,
                            resRange);

    // pull hists
    pull_prt[par] = new TH1F(Form("pull_prt_%s", par.c_str()),
                             Form("pull of %s", par.c_str()),
                             100,
                             -1 * pullRange,
                             pullRange);
    pull_flt[par] = new TH1F(Form("pull_flt_%s", par.c_str()),
                             Form("pull of %s", par.c_str()),
                             100,
                             -1 * pullRange,
                             pullRange);
    pull_smt[par] = new TH1F(Form("pull_smt_%s", par.c_str()),
                             Form("pull of %s", par.c_str()),
                             100,
                             -1 * pullRange,
                             pullRange);

    res_prt[par]->GetXaxis()->SetTitle(Form("r_{%s}", par.c_str()));
    res_prt[par]->GetYaxis()->SetTitle("Entries");
    res_flt[par]->GetXaxis()->SetTitle(Form("r_{%s}", par.c_str()));
    res_flt[par]->GetYaxis()->SetTitle("Entries");
    res_smt[par]->GetXaxis()->SetTitle(Form("r_{%s}", par.c_str()));
    res_smt[par]->GetYaxis()->SetTitle("Entries");

    pull_prt[par]->GetXaxis()->SetTitle(Form("pull_{%s}", par.c_str()));
    pull_prt[par]->GetYaxis()->SetTitle("Entries");
    pull_flt[par]->GetXaxis()->SetTitle(Form("pull_{%s}", par.c_str()));
    pull_flt[par]->GetYaxis()->SetTitle("Entries");
    pull_smt[par]->GetXaxis()->SetTitle(Form("pull_{%s}", par.c_str()));
    pull_smt[par]->GetYaxis()->SetTitle("Entries");

    // set style
    setHistStyle(res_prt[par], 1);
    setHistStyle(res_flt[par], 2);
    setHistStyle(res_smt[par], 4);

    setHistStyle(pull_prt[par], 1);
    setHistStyle(pull_flt[par], 2);
    setHistStyle(pull_smt[par], 4);
  }

  std::vector<float>* LOC0_prt
      = new std::vector<float>;  ///< predicted parameter local x
  std::vector<float>* LOC1_prt
      = new std::vector<float>;  ///< predicted parameter local y
  std::vector<float>* PHI_prt
      = new std::vector<float>;  ///< predicted parameter phi
  std::vector<float>* THETA_prt
      = new std::vector<float>;  ///< predicted parameter theta
  std::vector<float>* QOP_prt
      = new std::vector<float>;  ///< predicted parameter q/p
  std::vector<float>* T_prt
      = new std::vector<float>;  ///< predicted parameter t 
  std::vector<float>* LOC0_flt
      = new std::vector<float>;  ///< filtered parameter local x
  std::vector<float>* LOC1_flt
      = new std::vector<float>;  ///< filtered parameter local y
  std::vector<float>* PHI_flt
      = new std::vector<float>;  ///< filtered parameter phi
  std::vector<float>* THETA_flt
      = new std::vector<float>;  ///< filtered parameter theta
  std::vector<float>* QOP_flt
      = new std::vector<float>;  ///< filtered parameter q/p
  std::vector<float>* T_flt
      = new std::vector<float>;  ///< filtered parameter t 
  std::vector<float>* LOC0_smt
      = new std::vector<float>;  ///< smoothed parameter local x
  std::vector<float>* LOC1_smt
      = new std::vector<float>;  ///< smoothed parameter local y
  std::vector<float>* PHI_smt
      = new std::vector<float>;  ///< smoothed parameter phi
  std::vector<float>* THETA_smt
      = new std::vector<float>;  ///< smoothed parameter theta
  std::vector<float>* QOP_smt
      = new std::vector<float>;  ///< smoothed parameter q/p
  std::vector<float>* T_smt
      = new std::vector<float>;  ///< smoothed parameter t 

  std::vector<float>* res_LOC0_prt
      = new std::vector<float>;  ///< residual of predicted parameter local x
  std::vector<float>* res_LOC1_prt
      = new std::vector<float>;  ///< residual of predicted parameter local y
  std::vector<float>* res_PHI_prt
      = new std::vector<float>;  ///< residual of predicted parameter phi
  std::vector<float>* res_THETA_prt
      = new std::vector<float>;  ///< residual of predicted parameter theta
  std::vector<float>* res_QOP_prt
      = new std::vector<float>;  ///< residual of predicted parameter q/p
  std::vector<float>* res_T_prt
      = new std::vector<float>;  ///< residual of predicted parameter t 
  std::vector<float>* res_LOC0_flt
      = new std::vector<float>;  ///< residual of filtered parameter local x
  std::vector<float>* res_LOC1_flt
      = new std::vector<float>;  ///< residual of filtered parameter local y
  std::vector<float>* res_PHI_flt
      = new std::vector<float>;  ///< residual of filtered parameter phi
  std::vector<float>* res_THETA_flt
      = new std::vector<float>;  ///< residual of filtered parameter theta
  std::vector<float>* res_QOP_flt
      = new std::vector<float>;  ///< residual of filtered parameter q/p
  std::vector<float>* res_T_flt
      = new std::vector<float>;  ///< residual of filtered parameter t 
  std::vector<float>* res_LOC0_smt
      = new std::vector<float>;  ///< residual of smoothed parameter local x
  std::vector<float>* res_LOC1_smt
      = new std::vector<float>;  ///< residual of smoothed parameter local y
  std::vector<float>* res_PHI_smt
      = new std::vector<float>;  ///< residual of smoothed parameter phi
  std::vector<float>* res_THETA_smt
      = new std::vector<float>;  ///< residual of smoothed parameter theta
  std::vector<float>* res_QOP_smt
      = new std::vector<float>;  ///< residual of smoothed parameter q/p
  std::vector<float>* res_T_smt
      = new std::vector<float>;  ///< residual of smoothed parameter t 

  std::vector<float>* pull_LOC0_prt
      = new std::vector<float>;  ///< pull of predicted parameter local x
  std::vector<float>* pull_LOC1_prt
      = new std::vector<float>;  ///< pull of predicted parameter local y
  std::vector<float>* pull_PHI_prt
      = new std::vector<float>;  ///< pull of predicted parameter phi
  std::vector<float>* pull_THETA_prt
      = new std::vector<float>;  ///< pull of predicted parameter theta
  std::vector<float>* pull_QOP_prt
      = new std::vector<float>;  ///< pull of predicted parameter q/p
  std::vector<float>* pull_T_prt
      = new std::vector<float>;  ///< pull of predicted parameter t 
  std::vector<float>* pull_LOC0_flt
      = new std::vector<float>;  ///< pull of filtered parameter local x
  std::vector<float>* pull_LOC1_flt
      = new std::vector<float>;  ///< pull of filtered parameter local y
  std::vector<float>* pull_PHI_flt
      = new std::vector<float>;  ///< pull of filtered parameter phi
  std::vector<float>* pull_THETA_flt
      = new std::vector<float>;  ///< pull of filtered parameter theta
  std::vector<float>* pull_QOP_flt
      = new std::vector<float>;  ///< pull of filtered parameter q/p
  std::vector<float>* pull_T_flt
      = new std::vector<float>;  ///< pull of filtered parameter t 
  std::vector<float>* pull_LOC0_smt
      = new std::vector<float>;  ///< pull of smoothed parameter local x
  std::vector<float>* pull_LOC1_smt
      = new std::vector<float>;  ///< pull of smoothed parameter local y
  std::vector<float>* pull_PHI_smt
      = new std::vector<float>;  ///< pull of smoothed parameter phi
  std::vector<float>* pull_THETA_smt
      = new std::vector<float>;  ///< pull of smoothed parameter theta
  std::vector<float>* pull_QOP_smt
      = new std::vector<float>;  ///< pull of smoothed parameter q/p
  std::vector<float>* pull_T_smt
      = new std::vector<float>;  ///< pull of smoothed parameter t 

  std::vector<int>* volume_id = new std::vector<int>;  ///< volume_id
  std::vector<int>* layer_id  = new std::vector<int>;  ///< layer_id
  std::vector<int>* module_id = new std::vector<int>;  ///< module_id

  std::vector<bool>* predicted = new std::vector<bool>;  ///< prediction status
  std::vector<bool>* filtered  = new std::vector<bool>;  ///< filtering status
  std::vector<bool>* smoothed  = new std::vector<bool>;  ///< smoothing status

  int nStates, nMeasurements;

  tree->SetBranchAddress("eLOC0_prt", &LOC0_prt);
  tree->SetBranchAddress("eLOC1_prt", &LOC1_prt);
  tree->SetBranchAddress("ePHI_prt", &PHI_prt);
  tree->SetBranchAddress("eTHETA_prt", &THETA_prt);
  tree->SetBranchAddress("eQOP_prt", &QOP_prt);
  tree->SetBranchAddress("eT_prt", &T_prt);
  tree->SetBranchAddress("eLOC0_flt", &LOC0_flt);
  tree->SetBranchAddress("eLOC1_flt", &LOC1_flt);
  tree->SetBranchAddress("ePHI_flt", &PHI_flt);
  tree->SetBranchAddress("eTHETA_flt", &THETA_flt);
  tree->SetBranchAddress("eQOP_flt", &QOP_flt);
  tree->SetBranchAddress("eT_flt", &T_flt);
  tree->SetBranchAddress("eLOC0_smt", &LOC0_smt);
  tree->SetBranchAddress("eLOC1_smt", &LOC1_smt);
  tree->SetBranchAddress("ePHI_smt", &PHI_smt);
  tree->SetBranchAddress("eTHETA_smt", &THETA_smt);
  tree->SetBranchAddress("eQOP_smt", &QOP_smt);
  tree->SetBranchAddress("eT_smt", &T_smt);

  tree->SetBranchAddress("res_eLOC0_prt", &res_LOC0_prt);
  tree->SetBranchAddress("res_eLOC1_prt", &res_LOC1_prt);
  tree->SetBranchAddress("res_ePHI_prt", &res_PHI_prt);
  tree->SetBranchAddress("res_eTHETA_prt", &res_THETA_prt);
  tree->SetBranchAddress("res_eQOP_prt", &res_QOP_prt);
  tree->SetBranchAddress("res_eT_prt", &res_T_prt);
  tree->SetBranchAddress("res_eLOC0_flt", &res_LOC0_flt);
  tree->SetBranchAddress("res_eLOC1_flt", &res_LOC1_flt);
  tree->SetBranchAddress("res_ePHI_flt", &res_PHI_flt);
  tree->SetBranchAddress("res_eTHETA_flt", &res_THETA_flt);
  tree->SetBranchAddress("res_eQOP_flt", &res_QOP_flt);
  tree->SetBranchAddress("res_eT_flt", &res_T_flt);
  tree->SetBranchAddress("res_eLOC0_smt", &res_LOC0_smt);
  tree->SetBranchAddress("res_eLOC1_smt", &res_LOC1_smt);
  tree->SetBranchAddress("res_ePHI_smt", &res_PHI_smt);
  tree->SetBranchAddress("res_eTHETA_smt", &res_THETA_smt);
  tree->SetBranchAddress("res_eQOP_smt", &res_QOP_smt);
  tree->SetBranchAddress("res_eT_smt", &res_T_smt);

  tree->SetBranchAddress("pull_eLOC0_prt", &pull_LOC0_prt);
  tree->SetBranchAddress("pull_eLOC1_prt", &pull_LOC1_prt);
  tree->SetBranchAddress("pull_ePHI_prt", &pull_PHI_prt);
  tree->SetBranchAddress("pull_eTHETA_prt", &pull_THETA_prt);
  tree->SetBranchAddress("pull_eQOP_prt", &pull_QOP_prt);
  tree->SetBranchAddress("pull_eT_prt", &pull_T_prt);
  tree->SetBranchAddress("pull_eLOC0_flt", &pull_LOC0_flt);
  tree->SetBranchAddress("pull_eLOC1_flt", &pull_LOC1_flt);
  tree->SetBranchAddress("pull_ePHI_flt", &pull_PHI_flt);
  tree->SetBranchAddress("pull_eTHETA_flt", &pull_THETA_flt);
  tree->SetBranchAddress("pull_eQOP_flt", &pull_QOP_flt);
  tree->SetBranchAddress("pull_eT_flt", &pull_T_flt);
  tree->SetBranchAddress("pull_eLOC0_smt", &pull_LOC0_smt);
  tree->SetBranchAddress("pull_eLOC1_smt", &pull_LOC1_smt);
  tree->SetBranchAddress("pull_ePHI_smt", &pull_PHI_smt);
  tree->SetBranchAddress("pull_eTHETA_smt", &pull_THETA_smt);
  tree->SetBranchAddress("pull_eQOP_smt", &pull_QOP_smt);
  tree->SetBranchAddress("pull_eT_smt", &pull_T_smt);

  tree->SetBranchAddress("nStates", &nStates);
  tree->SetBranchAddress("nMeasurements", &nMeasurements);
  tree->SetBranchAddress("volume_id", &volume_id);
  tree->SetBranchAddress("layer_id", &layer_id);
  tree->SetBranchAddress("module_id", &module_id);
  tree->SetBranchAddress("predicted", &predicted);
  tree->SetBranchAddress("filtered", &filtered);
  tree->SetBranchAddress("smoothed", &smoothed);

  Int_t entries = tree->GetEntries();
  for (int j = 0; j < entries; j++) {
    tree->GetEvent(j);
  
    for (int i = 0; i < nMeasurements; i++) {
      if (predicted->at(i)) {
        res_prt[paramNames[0]]->Fill(res_LOC0_prt->at(i), 1);
        res_prt[paramNames[1]]->Fill(res_LOC1_prt->at(i), 1);
        res_prt[paramNames[2]]->Fill(res_PHI_prt->at(i), 1);
        res_prt[paramNames[3]]->Fill(res_THETA_prt->at(i), 1);
        res_prt[paramNames[4]]->Fill(res_QOP_prt->at(i), 1);
        res_prt[paramNames[5]]->Fill(res_T_prt->at(i), 1);
        pull_prt[paramNames[0]]->Fill(pull_LOC0_prt->at(i), 1);
        pull_prt[paramNames[1]]->Fill(pull_LOC1_prt->at(i), 1);
        pull_prt[paramNames[2]]->Fill(pull_PHI_prt->at(i), 1);
        pull_prt[paramNames[3]]->Fill(pull_THETA_prt->at(i), 1);
        pull_prt[paramNames[4]]->Fill(pull_QOP_prt->at(i), 1);
        pull_prt[paramNames[5]]->Fill(pull_T_prt->at(i), 1);
      }
      if (filtered->at(i)) {
        res_flt[paramNames[0]]->Fill(res_LOC0_flt->at(i), 1);
        res_flt[paramNames[1]]->Fill(res_LOC1_flt->at(i), 1);
        res_flt[paramNames[2]]->Fill(res_PHI_flt->at(i), 1);
        res_flt[paramNames[3]]->Fill(res_THETA_flt->at(i), 1);
        res_flt[paramNames[4]]->Fill(res_QOP_flt->at(i), 1);
        res_flt[paramNames[5]]->Fill(res_T_flt->at(i), 1);
        pull_flt[paramNames[0]]->Fill(pull_LOC0_flt->at(i), 1);
        pull_flt[paramNames[1]]->Fill(pull_LOC1_flt->at(i), 1);
        pull_flt[paramNames[2]]->Fill(pull_PHI_flt->at(i), 1);
        pull_flt[paramNames[3]]->Fill(pull_THETA_flt->at(i), 1);
        pull_flt[paramNames[4]]->Fill(pull_QOP_flt->at(i), 1);
        pull_flt[paramNames[5]]->Fill(pull_T_flt->at(i), 1);
      }
      if (smoothed->at(i)) {
        res_smt[paramNames[0]]->Fill(res_LOC0_smt->at(i), 1);
        res_smt[paramNames[1]]->Fill(res_LOC1_smt->at(i), 1);
        res_smt[paramNames[2]]->Fill(res_PHI_smt->at(i), 1);
        res_smt[paramNames[3]]->Fill(res_THETA_smt->at(i), 1);
        res_smt[paramNames[4]]->Fill(res_QOP_smt->at(i), 1);
        res_smt[paramNames[5]]->Fill(res_T_smt->at(i), 1);
        pull_smt[paramNames[0]]->Fill(pull_LOC0_smt->at(i), 1);
        pull_smt[paramNames[1]]->Fill(pull_LOC1_smt->at(i), 1);
        pull_smt[paramNames[2]]->Fill(pull_PHI_smt->at(i), 1);
        pull_smt[paramNames[3]]->Fill(pull_THETA_smt->at(i), 1);
        pull_smt[paramNames[4]]->Fill(pull_QOP_smt->at(i), 1);
        pull_smt[paramNames[5]]->Fill(pull_T_smt->at(i), 1);
      }
    }
  }

  // plotting residual
  TCanvas* c1 = new TCanvas("c1", "c1", 1200, 800);
  c1->Divide(3, 2);
  for (size_t ipar = 0; ipar < paramNames.size(); ipar++) {
    c1->cd(ipar + 1);
    res_smt[paramNames.at(ipar)]->Draw("");
    res_prt[paramNames.at(ipar)]->Draw("same");
    res_flt[paramNames.at(ipar)]->Draw("same");

    int binmax     = res_smt[paramNames.at(ipar)]->GetMaximumBin();
    int bincontent = res_smt[paramNames.at(ipar)]->GetBinContent(binmax);

    res_smt[paramNames.at(ipar)]->GetYaxis()->SetRangeUser(0, bincontent * 1.2);
    TLegend* legend = new TLegend(0.7, 0.7, 0.9, 0.9);
    legend->AddEntry(res_prt[paramNames.at(ipar)], "prediction", "lp");
    legend->AddEntry(res_flt[paramNames.at(ipar)], "filtering", "lp");
    legend->AddEntry(res_smt[paramNames.at(ipar)], "smoothing", "lp");
    legend->SetBorderSize(0);
    legend->SetFillStyle(0);
    legend->SetTextFont(42);
    legend->Draw();
  }

  // plotting pull
  TCanvas* c2 = new TCanvas("c2", "c2", 1200, 800);
  c2->Divide(3, 2);
  for (size_t ipar = 0; ipar < paramNames.size(); ipar++) {
    c2->cd(ipar + 1);
    pull_smt[paramNames.at(ipar)]->Draw("");
    pull_prt[paramNames.at(ipar)]->Draw("same");
    pull_flt[paramNames.at(ipar)]->Draw("same");

    int binmax     = pull_smt[paramNames.at(ipar)]->GetMaximumBin();
    int bincontent = pull_smt[paramNames.at(ipar)]->GetBinContent(binmax);

    pull_smt[paramNames.at(ipar)]->GetYaxis()->SetRangeUser(0,
                                                            bincontent * 1.2);
    TLegend* legend = new TLegend(0.7, 0.7, 0.9, 0.9);
    legend->AddEntry(pull_prt[paramNames.at(ipar)], "prediction", "lp");
    legend->AddEntry(pull_flt[paramNames.at(ipar)], "filtering", "lp");
    legend->AddEntry(pull_smt[paramNames.at(ipar)], "smoothing", "lp");
    legend->SetBorderSize(0);
    legend->SetFillStyle(0);
    legend->SetTextFont(42);
    legend->Draw();
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
  hist->SetLineColor(color);
  hist->SetMarkerColor(color);
}

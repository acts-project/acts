// This file is part of the Acts project.
//
// Copyright (C) 2019-2021 CERN for the benefit of the Acts project
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
#include <TF1.h>
#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TLegend.h>
#include <TMath.h>
#include <TProfile2D.h>
#include <TStyle.h>
#include <TTree.h>

#include "CommonUtils.h"
#include "TreeReader.h"

using namespace ROOT;

// (loc1, phi, theta, q/p, t) at all track states from root file produced by the
// RootTrajectoryStatesWriter
//
void boundParamResolution(const std::string& inFile,
                          const std::string& treeName,
                          const std::string& outFile,
                          const std::string& saveAs = "") {
  gStyle->SetOptFit(0000);
  gStyle->SetOptStat(0000);
  gStyle->SetPadLeftMargin(0.20);
  gStyle->SetPadRightMargin(0.1);
  gStyle->SetPadTopMargin(0.1);
  gStyle->SetPadBottomMargin(0.15);

  // Pull range
  double pullRange = 7;

  // Section 0: file handling ---------------------------------------------
  //
  // Open root file written by RootTrajectoryWriter
  // Create output root file
  std::cout << "Opening file: " << inFile << std::endl;
  TFile* file = TFile::Open(inFile.c_str(), "read");
  std::cout << "Reading tree: " << treeName << std::endl;
  TTree* tree = (TTree*)file->Get(treeName.c_str());
  TFile* out = TFile::Open(outFile.c_str(), "recreate");
  out->cd();

  // Section 1: geometry parsing -------------------------------------------
  //
  // Gathers accessible volume_id and layer_id values
  // Draw the volume_id : layer_id correlation matrix
  std::map<int, std::vector<int>> volLayIds;
  TCanvas* geometryCanvas =
      new TCanvas("volLayCanvas", "Volume Layer Matrix", 10, 10, 450, 600);
  geometryCanvas->Divide(1, 2);
  // Volume/Layer Id
  geometryCanvas->cd(1);
  tree->Draw("layer_id:volume_id>>volID_layID(100,0.5,100.5,100,0.5,100.5)", "",
             "colz");
  auto h2_volID_layID = dynamic_cast<TH2F*>(out->Get("volID_layID"));
  setHistStyle(h2_volID_layID, 2);
  h2_volID_layID->Write();
  // Geometry size
  geometryCanvas->cd(2);
  tree->Draw("sqrt(g_x_smt*g_x_smt+g_y_smt*g_y_smt):g_z_smt>>geo_dim");
  auto h2_geo_dim = dynamic_cast<TH2F*>(out->Get("geo_dim"));
  setHistStyle(h2_geo_dim, 2);
  float detectorZ = 1.15 * h2_geo_dim->GetXaxis()->GetXmax();
  float detectorR = 1.15 * h2_geo_dim->GetYaxis()->GetXmax();
  h2_geo_dim->Write();

  // Save the plots on screen
  if (not saveAs.empty()) {
    geometryCanvas->SaveAs((std::string("all_vol_lay_ids.") + saveAs).c_str());
  }

  // Now determine the valid volume / layer bins
  int volBins = h2_volID_layID->GetXaxis()->GetNbins();
  int layBins = h2_volID_layID->GetYaxis()->GetNbins();
  // Add the overall plots
  volLayIds[-1] = {-1};
  // Go through volumes and add plots per volume
  for (int volID = 1; volID <= volBins; ++volID) {
    for (int layID = 1; layID <= layBins; ++layID) {
      if (h2_volID_layID->GetBinContent(volID, layID) != 0.) {
        if (volLayIds.find(volID) == volLayIds.end()) {
          // First occurance of this layer, add -1 for volume plots
          volLayIds[volID] = {-1, layID};
        } else {
          volLayIds[volID].push_back(layID);
        }
      }
    }
  }

  // Section 2: Branch assignment
  //
  // Helper for assigning branches to the input tree
  TrackStatesReader tsReader(tree, false);


  // Section 3: Histogram booking
  TCanvas* rangeCanvas =
      new TCanvas("rangeCanvas", "Range Estimation", 10, 10, 900, 600);
  rangeCanvas->Divide(3, 2);
  std::vector<std::string> res_smts = {"res_eLOC0_smt", "res_eLOC1_smt",
                                       "res_ePHI_smt",  "res_eTHETA_smt",
                                       "res_eQOP_smt",  "res_eT_smt"};

  /// Helper method to make a volLayId string for identification & cut
  ///
  /// ensures a unique identification
  auto volLayIdCut = [](int vol, int lay) -> std::array<std::string, 2> {
    if (vol < 0 and lay < 0) {
      return {std::string("all_"), ""};
    }

    std::string vlstr = "vol" + std::to_string(vol);
    std::string vlcut = "volume_id == " + std::to_string(vol);
    if (lay > 0) {
      vlstr += (std::string("_lay") + std::to_string(lay));
      vlcut += (std::string(" && layer_id == ") + std::to_string(lay));
    }
    vlstr += std::string("_");
    return {vlstr, vlcut};
  };

  /// Helper method to estimate the ranges
  ///
  /// Takes identification & cut from the first one
  auto histRanges =
      [&](const std::array<std::string, 2>& vlIdCut) -> std::array<float, 6> {
    std::array<float, 6> ranges = {0., 0., 0., 0., 0., 0.};
    for (unsigned int ir = 0; ir < 6; ++ir) {
      rangeCanvas->cd(ir + 1);
      std::string drawString = res_smts[ir] + ">>";
      std::string histString =
          vlIdCut[0] + std::string("range_") + res_smts[ir];
      tree->Draw((drawString + histString).c_str(), vlIdCut[1].c_str());
      auto h1_range = dynamic_cast<TH1F*>(out->Get(histString.c_str()));
      h1_range->Write();
      float range = pullRange * h1_range->GetRMS();
      ranges[ir] = range;
    }
    if (not saveAs.empty()) {
      rangeCanvas->SaveAs(
          (vlIdCut[0] + std::string("res_ranges.") + saveAs).c_str());
    }
    return ranges;
  };

  // Track parameter name
  std::vector<std::string> paramNames = {"loc0",   "loc1", "#phi",
                                         "#theta", "q/p",  "t"};

  // Book histograms (with adapted range):
  //
  // Global profile histograms : residuals/pulls
  std::array<TProfile2D*, 6> p2d_res_zr_prt;
  std::array<TProfile2D*, 6> p2d_res_zr_flt;
  std::array<TProfile2D*, 6> p2d_res_zr_smt;
  std::array<TProfile2D*, 6> p2d_pull_zr_prt;
  std::array<TProfile2D*, 6> p2d_pull_zr_flt;
  std::array<TProfile2D*, 6> p2d_pull_zr_smt;

  for (unsigned int ipar = 0; ipar < paramNames.size(); ++ipar) {
    const auto& par = paramNames[ipar];
    p2d_res_zr_prt[ipar] =
        new TProfile2D(Form("all_prof_res_prt_%s", par.c_str()),
                       Form("residual profile of %s", par.c_str()), 100,
                       -detectorZ, detectorZ, 50, 0., detectorR);
    p2d_res_zr_flt[ipar] =
        new TProfile2D(Form("all_prof_res_flt_%s", par.c_str()),
                       Form("residual profile of %s", par.c_str()), 100,
                       -detectorZ, detectorZ, 50, 0., detectorR);
    p2d_res_zr_smt[ipar] =
        new TProfile2D(Form("all_prof_res_smt_%s", par.c_str()),
                       Form("residual profile of %s", par.c_str()), 100,
                       -detectorZ, detectorZ, 50, 0., detectorR);

    p2d_pull_zr_prt[ipar] =
        new TProfile2D(Form("all_prof_pull_prt_%s", par.c_str()),
                       Form("pull profile of %s", par.c_str()), 100,
                       -detectorZ, detectorZ, 50, 0., detectorR);
    p2d_pull_zr_flt[ipar] =
        new TProfile2D(Form("all_prof_pull_flt_%s", par.c_str()),
                       Form("pull profile of %s", par.c_str()), 100,
                       -detectorZ, detectorZ, 50, 0., detectorR);
    p2d_pull_zr_smt[ipar] =
        new TProfile2D(Form("all_prof_pull_smt_%s", par.c_str()),
                       Form("pull profile of %s", par.c_str()), 100,
                       -detectorZ, detectorZ, 50, 0., detectorR);

    p2d_res_zr_prt[ipar]->SetErrorOption("s");
    p2d_res_zr_flt[ipar]->SetErrorOption("s");
    p2d_res_zr_smt[ipar]->SetErrorOption("s");
    p2d_pull_zr_prt[ipar]->SetErrorOption("s");
    p2d_pull_zr_flt[ipar]->SetErrorOption("s");
    p2d_pull_zr_smt[ipar]->SetErrorOption("s");

    p2d_res_zr_prt[ipar]->GetXaxis()->SetTitle("z [mm]");
    p2d_res_zr_prt[ipar]->GetYaxis()->SetTitle("R [mm]");
    p2d_res_zr_prt[ipar]->GetZaxis()->SetTitle(Form("r_{%s}", par.c_str()));
    p2d_res_zr_flt[ipar]->GetXaxis()->SetTitle("z [mm]");
    p2d_res_zr_flt[ipar]->GetYaxis()->SetTitle("R [mm]");
    p2d_res_zr_flt[ipar]->GetZaxis()->SetTitle(Form("r_{%s}", par.c_str()));
    p2d_res_zr_smt[ipar]->GetXaxis()->SetTitle("z [mm]");
    p2d_res_zr_smt[ipar]->GetYaxis()->SetTitle("R [mm]");
    p2d_res_zr_smt[ipar]->GetZaxis()->SetTitle(Form("r_{%s}", par.c_str()));

    p2d_pull_zr_prt[ipar]->GetXaxis()->SetTitle("z [mm]");
    p2d_pull_zr_prt[ipar]->GetYaxis()->SetTitle("R [mm]");
    p2d_pull_zr_prt[ipar]->GetZaxis()->SetTitle(Form("pull_{%s}", par.c_str()));
    p2d_pull_zr_flt[ipar]->GetXaxis()->SetTitle("z [mm]");
    p2d_pull_zr_flt[ipar]->GetYaxis()->SetTitle("R [mm]");
    p2d_pull_zr_prt[ipar]->GetZaxis()->SetTitle(Form("pull_{%s}", par.c_str()));
    p2d_pull_zr_smt[ipar]->GetXaxis()->SetTitle("z [mm]");
    p2d_pull_zr_smt[ipar]->GetYaxis()->SetTitle("R [mm]");
    p2d_pull_zr_prt[ipar]->GetZaxis()->SetTitle(Form("pull_{%s}", par.c_str()));

    // set style
    setHistStyle(p2d_res_zr_prt[ipar], 1);
    setHistStyle(p2d_res_zr_flt[ipar], 2);
    setHistStyle(p2d_res_zr_smt[ipar], 4);

    setHistStyle(p2d_pull_zr_prt[ipar], 1);
    setHistStyle(p2d_pull_zr_flt[ipar], 2);
    setHistStyle(p2d_pull_zr_smt[ipar], 4);
  }

  // Resiudal / Pull histograms
  std::map<std::string, TH1F*> res_prt;
  std::map<std::string, TH1F*> res_flt;
  std::map<std::string, TH1F*> res_smt;
  std::map<std::string, TH1F*> pull_prt;
  std::map<std::string, TH1F*> pull_flt;
  std::map<std::string, TH1F*> pull_smt;

  // - per layer (identified by vol, lay)
  // - per volume (identified by vol, -1)
  // - overall (identified by lay)
  for (auto [vol, layers] : volLayIds) {
    for (auto lay : layers) {
      // Estimate the ranges from smoothed
      auto vlIdCut = volLayIdCut(vol, lay);
      auto ranges = histRanges(vlIdCut);

      // Residual range
      std::map<std::string, double> paramResidualRange = {
          {"loc0", ranges[0]},   {"loc1", ranges[1]}, {"#phi", ranges[2]},
          {"#theta", ranges[3]}, {"q/p", ranges[4]},  {"t", ranges[5]}};

      // Create the hists and set up for them
      for (const auto& [par, resRange] : paramResidualRange) {
        // histogram creation
        std::string id_par = vlIdCut[0] + par;
        // residual hists
        res_prt[id_par] = new TH1F(
            Form((vlIdCut[0] + std::string("res_prt_%s")).c_str(), par.c_str()),
            Form("residual of %s", par.c_str()), 100, -1 * resRange, resRange);
        res_flt[id_par] = new TH1F(
            Form((vlIdCut[0] + std::string("res_flt_%s")).c_str(), par.c_str()),
            Form("residual of %s", par.c_str()), 100, -1 * resRange, resRange);
        res_smt[id_par] = new TH1F(
            Form((vlIdCut[0] + std::string("res_smt_%s")).c_str(), par.c_str()),
            Form("residual of %s", par.c_str()), 100, -1 * resRange, resRange);

        // pull hists
        pull_prt[id_par] = new TH1F(
            Form((vlIdCut[0] + std::string("pull_prt_%s")).c_str(),
                 par.c_str()),
            Form("pull of %s", par.c_str()), 100, -1 * pullRange, pullRange);
        pull_flt[id_par] = new TH1F(
            Form((vlIdCut[0] + std::string("pull_flt_%s")).c_str(),
                 par.c_str()),
            Form("pull of %s", par.c_str()), 100, -1 * pullRange, pullRange);
        pull_smt[id_par] = new TH1F(
            Form((vlIdCut[0] + std::string("pull_smt_%s")).c_str(),
                 par.c_str()),
            Form("pull of %s", par.c_str()), 100, -1 * pullRange, pullRange);

        res_prt[id_par]->GetXaxis()->SetTitle(Form("r_{%s}", par.c_str()));
        res_prt[id_par]->GetYaxis()->SetTitle("Entries");
        res_flt[id_par]->GetXaxis()->SetTitle(Form("r_{%s}", par.c_str()));
        res_flt[id_par]->GetYaxis()->SetTitle("Entries");
        res_smt[id_par]->GetXaxis()->SetTitle(Form("r_{%s}", par.c_str()));
        res_smt[id_par]->GetYaxis()->SetTitle("Entries");

        pull_prt[id_par]->GetXaxis()->SetTitle(Form("pull_{%s}", par.c_str()));
        pull_prt[id_par]->GetYaxis()->SetTitle("Entries");
        pull_flt[id_par]->GetXaxis()->SetTitle(Form("pull_{%s}", par.c_str()));
        pull_flt[id_par]->GetYaxis()->SetTitle("Entries");
        pull_smt[id_par]->GetXaxis()->SetTitle(Form("pull_{%s}", par.c_str()));
        pull_smt[id_par]->GetYaxis()->SetTitle("Entries");

        // set style
        setHistStyle(res_prt[id_par], 1);
        setHistStyle(res_flt[id_par], 2);
        setHistStyle(res_smt[id_par], 4);

        setHistStyle(pull_prt[id_par], 1);
        setHistStyle(pull_flt[id_par], 2);
        setHistStyle(pull_smt[id_par], 4);
      }
    }
  }

  // Section 4: Histogram filling
  //
  // - Running through the entries and filling the histograms
  int entries = tree->GetEntries();
  for (int j = 0; j < entries; j++) {
    tsReader.getEntry(j);

    for (unsigned int i = 0; i < tsReader.nMeasurements; i++) {
      // global profile filling
      if (tsReader.predicted->at(i)) {
        float x_prt = tsReader.g_x_prt->at(i);
        float y_prt = tsReader.g_y_prt->at(i);
        float r_prt = std::sqrt(x_prt * x_prt + y_prt * y_prt);
        float z_prt = tsReader.g_z_prt->at(i);
        p2d_res_zr_prt[0]->Fill(z_prt, r_prt, tsReader.res_LOC0_prt->at(i));
        p2d_res_zr_prt[1]->Fill(z_prt, r_prt, tsReader.res_LOC1_prt->at(i));
        p2d_res_zr_prt[2]->Fill(z_prt, r_prt, tsReader.res_PHI_prt->at(i));
        p2d_res_zr_prt[3]->Fill(z_prt, r_prt, tsReader.res_THETA_prt->at(i));
        p2d_res_zr_prt[4]->Fill(z_prt, r_prt, tsReader.res_QOP_prt->at(i));
        p2d_res_zr_prt[5]->Fill(z_prt, r_prt, tsReader.res_T_prt->at(i));
        p2d_pull_zr_prt[0]->Fill(z_prt, r_prt, tsReader.pull_LOC0_prt->at(i));
        p2d_pull_zr_prt[1]->Fill(z_prt, r_prt, tsReader.pull_LOC1_prt->at(i));
        p2d_pull_zr_prt[2]->Fill(z_prt, r_prt, tsReader.pull_PHI_prt->at(i));
        p2d_pull_zr_prt[3]->Fill(z_prt, r_prt, tsReader.pull_THETA_prt->at(i));
        p2d_pull_zr_prt[4]->Fill(z_prt, r_prt, tsReader.pull_QOP_prt->at(i));
        p2d_pull_zr_prt[5]->Fill(z_prt, r_prt, tsReader.pull_T_prt->at(i));
      }
      if (tsReader.filtered->at(i)) {
        float x_flt = tsReader.g_x_flt->at(i);
        float y_flt = tsReader.g_y_flt->at(i);
        float r_flt = std::sqrt(x_flt * x_flt + y_flt * y_flt);
        float z_flt = tsReader.g_z_flt->at(i);
        p2d_res_zr_flt[0]->Fill(z_flt, r_flt, tsReader.res_LOC0_flt->at(i));
        p2d_res_zr_flt[1]->Fill(z_flt, r_flt, tsReader.res_LOC1_flt->at(i));
        p2d_res_zr_flt[2]->Fill(z_flt, r_flt, tsReader.res_PHI_flt->at(i));
        p2d_res_zr_flt[3]->Fill(z_flt, r_flt, tsReader.res_THETA_flt->at(i));
        p2d_res_zr_flt[4]->Fill(z_flt, r_flt, tsReader.res_QOP_flt->at(i));
        p2d_res_zr_flt[5]->Fill(z_flt, r_flt, tsReader.res_T_flt->at(i));
        p2d_pull_zr_flt[0]->Fill(z_flt, r_flt, tsReader.pull_LOC0_flt->at(i));
        p2d_pull_zr_flt[1]->Fill(z_flt, r_flt, tsReader.pull_LOC1_flt->at(i));
        p2d_pull_zr_flt[2]->Fill(z_flt, r_flt, tsReader.pull_PHI_flt->at(i));
        p2d_pull_zr_flt[3]->Fill(z_flt, r_flt, tsReader.pull_THETA_flt->at(i));
        p2d_pull_zr_flt[4]->Fill(z_flt, r_flt, tsReader.pull_QOP_flt->at(i));
        p2d_pull_zr_flt[5]->Fill(z_flt, r_flt, tsReader.pull_T_flt->at(i));
      }
      if (tsReader.smoothed->at(i)) {
        float x_smt = tsReader.g_x_smt->at(i);
        float y_smt = tsReader.g_y_smt->at(i);
        float r_smt = std::sqrt(x_smt * x_smt + y_smt * y_smt);
        float z_smt = tsReader.g_z_smt->at(i);
        p2d_res_zr_smt[0]->Fill(z_smt, r_smt, tsReader.res_LOC0_smt->at(i));
        p2d_res_zr_smt[1]->Fill(z_smt, r_smt, tsReader.res_LOC1_smt->at(i));
        p2d_res_zr_smt[2]->Fill(z_smt, r_smt, tsReader.res_PHI_smt->at(i));
        p2d_res_zr_smt[3]->Fill(z_smt, r_smt, tsReader.res_THETA_smt->at(i));
        p2d_res_zr_smt[4]->Fill(z_smt, r_smt, tsReader.res_QOP_smt->at(i));
        p2d_res_zr_smt[5]->Fill(z_smt, r_smt, tsReader.res_T_smt->at(i));
        p2d_pull_zr_smt[0]->Fill(z_smt, r_smt, tsReader.pull_LOC0_smt->at(i));
        p2d_pull_zr_smt[1]->Fill(z_smt, r_smt, tsReader.pull_LOC1_smt->at(i));
        p2d_pull_zr_smt[2]->Fill(z_smt, r_smt, tsReader.pull_PHI_smt->at(i));
        p2d_pull_zr_smt[3]->Fill(z_smt, r_smt, tsReader.pull_THETA_smt->at(i));
        p2d_pull_zr_smt[4]->Fill(z_smt, r_smt, tsReader.pull_QOP_smt->at(i));
        p2d_pull_zr_smt[5]->Fill(z_smt, r_smt, tsReader.pull_T_smt->at(i));
      }

      int vol = tsReader.volume_id->at(i);
      int lay = tsReader.layer_id->at(i);

      /// Always fill (-1,-1), (vol, -1), (vol, lay)
      std::vector<std::array<int, 2>> fillIds = {
          {-1, -1}, {vol, -1}, {vol, lay}};

      for (const auto& fid : fillIds) {
        auto vlID = volLayIdCut(fid[0], fid[1])[0];
        // Fill predicated parameters
        if (tsReader.predicted->at(i)) {
          res_prt[vlID + paramNames[0]]->Fill(tsReader.res_LOC0_prt->at(i), 1);
          res_prt[vlID + paramNames[1]]->Fill(tsReader.res_LOC1_prt->at(i), 1);
          res_prt[vlID + paramNames[2]]->Fill(tsReader.res_PHI_prt->at(i), 1);
          res_prt[vlID + paramNames[3]]->Fill(tsReader.res_THETA_prt->at(i), 1);
          res_prt[vlID + paramNames[4]]->Fill(tsReader.res_QOP_prt->at(i), 1);
          res_prt[vlID + paramNames[5]]->Fill(tsReader.res_T_prt->at(i), 1);
          pull_prt[vlID + paramNames[0]]->Fill(tsReader.pull_LOC0_prt->at(i), 1);
          pull_prt[vlID + paramNames[1]]->Fill(tsReader.pull_LOC1_prt->at(i), 1);
          pull_prt[vlID + paramNames[2]]->Fill(tsReader.pull_PHI_prt->at(i), 1);
          pull_prt[vlID + paramNames[3]]->Fill(tsReader.pull_THETA_prt->at(i), 1);
          pull_prt[vlID + paramNames[4]]->Fill(tsReader.pull_QOP_prt->at(i), 1);
          pull_prt[vlID + paramNames[5]]->Fill(tsReader.pull_T_prt->at(i), 1);
        }
        // Fill filtered parameters
        if (tsReader.filtered->at(i)) {
          res_flt[vlID + paramNames[0]]->Fill(tsReader.res_LOC0_flt->at(i), 1);
          res_flt[vlID + paramNames[1]]->Fill(tsReader.res_LOC1_flt->at(i), 1);
          res_flt[vlID + paramNames[2]]->Fill(tsReader.res_PHI_flt->at(i), 1);
          res_flt[vlID + paramNames[3]]->Fill(tsReader.res_THETA_flt->at(i), 1);
          res_flt[vlID + paramNames[4]]->Fill(tsReader.res_QOP_flt->at(i), 1);
          res_flt[vlID + paramNames[5]]->Fill(tsReader.res_T_flt->at(i), 1);
          pull_flt[vlID + paramNames[0]]->Fill(tsReader.pull_LOC0_flt->at(i), 1);
          pull_flt[vlID + paramNames[1]]->Fill(tsReader.pull_LOC1_flt->at(i), 1);
          pull_flt[vlID + paramNames[2]]->Fill(tsReader.pull_PHI_flt->at(i), 1);
          pull_flt[vlID + paramNames[3]]->Fill(tsReader.pull_THETA_flt->at(i), 1);
          pull_flt[vlID + paramNames[4]]->Fill(tsReader.pull_QOP_flt->at(i), 1);
          pull_flt[vlID + paramNames[5]]->Fill(tsReader.pull_T_flt->at(i), 1);
        }
        // Fill smoothed parameters
        if (tsReader.smoothed->at(i)) {
          res_smt[vlID + paramNames[0]]->Fill(tsReader.res_LOC0_smt->at(i), 1);
          res_smt[vlID + paramNames[1]]->Fill(tsReader.res_LOC1_smt->at(i), 1);
          res_smt[vlID + paramNames[2]]->Fill(tsReader.res_PHI_smt->at(i), 1);
          res_smt[vlID + paramNames[3]]->Fill(tsReader.res_THETA_smt->at(i), 1);
          res_smt[vlID + paramNames[4]]->Fill(tsReader.res_QOP_smt->at(i), 1);
          res_smt[vlID + paramNames[5]]->Fill(tsReader.res_T_smt->at(i), 1);
          pull_smt[vlID + paramNames[0]]->Fill(tsReader.pull_LOC0_smt->at(i), 1);
          pull_smt[vlID + paramNames[1]]->Fill(tsReader.pull_LOC1_smt->at(i), 1);
          pull_smt[vlID + paramNames[2]]->Fill(tsReader.pull_PHI_smt->at(i), 1);
          pull_smt[vlID + paramNames[3]]->Fill(tsReader.pull_THETA_smt->at(i), 1);
          pull_smt[vlID + paramNames[4]]->Fill(tsReader.pull_QOP_smt->at(i), 1);
          pull_smt[vlID + paramNames[5]]->Fill(tsReader.pull_T_smt->at(i), 1);
        }
      }
    }
  }

  // Section 5: Histogram plotting

  // Plotting global profiles
  TCanvas* respull_mean_prf =
      new TCanvas("respull_mean_prf",
                  "Residual/Pull Distributions - mean profiles", 1800, 800);
  respull_mean_prf->Divide(3, 2);

  TCanvas* respull_var_prf =
      new TCanvas("respull_var_prf",
                  "Residual/Pull Distributions - variance profiles", 1800, 800);
  respull_var_prf->Divide(3, 2);

  auto plotProfiles = [&](std::array<TProfile2D*, 6>& profiles,
                          const std::string& res_pull,
                          const std::string& type) -> void {
    // Mean
    for (size_t ipar = 0; ipar < paramNames.size(); ++ipar) {
      respull_mean_prf->cd(ipar + 1);

      if (res_pull == "pull") {
        profiles[ipar]->GetZaxis()->SetRangeUser(-1. * pullRange, pullRange);
      }
      adaptColorPalette(profiles[ipar], -1. * pullRange, pullRange, 0., 0.25, 104);
      profiles[ipar]->Draw("colz");
      profiles[ipar]->Write();
    }
    // Save the canvas: mean
    if (not saveAs.empty()) {
      respull_mean_prf->SaveAs((std::string("all_") + res_pull +
                                std::string("_mean_prf_") + type +
                                std::string(".") + saveAs)
                                   .c_str());
    }

    // Variance
    for (size_t ipar = 0; ipar < paramNames.size(); ++ipar) {
      respull_var_prf->cd(ipar + 1);
      auto zAxis = profiles[ipar]->GetXaxis();
      auto rAxis = profiles[ipar]->GetYaxis();
      int binsZ = zAxis->GetNbins();
      int binsR = rAxis->GetNbins();
      std::string hist_name = "all_";
      hist_name += res_pull;
      hist_name += "_";
      hist_name += paramNames[ipar];
      hist_name += "_";
      hist_name += type;
      TH2F* var =
          new TH2F(hist_name.c_str(),
                   (res_pull + std::string(" ") + paramNames[ipar]).c_str(),
                   binsZ, zAxis->GetXmin(), zAxis->GetXmax(), binsR,
                   rAxis->GetXmin(), rAxis->GetXmax());
      for (int iz = 1; iz <= binsZ; ++iz) {
        for (int ir = 1; ir <= binsR; ++ir) {
          if (profiles[ipar]->GetBinContent(iz, ir) != 0.) {
            var->SetBinContent(iz, ir, profiles[ipar]->GetBinError(iz, ir));
          }
        }
      }
      if (res_pull == "pull") {
        var->GetZaxis()->SetRangeUser(0., pullRange);
      }
      adaptColorPalette(var, 0., pullRange, 1., 0.5, 104);
      var->Draw("colz");
      var->Write();
    }
    // Save the canvas: pulls
    if (not saveAs.empty()) {
      respull_var_prf->SaveAs((std::string("all_") + res_pull +
                               std::string("_var_prf_") + type +
                               std::string(".") + saveAs)
                                  .c_str());
    }
  };

  // Plotting profiles: res/pulls
  plotProfiles(p2d_res_zr_prt, "res", "prt");
  plotProfiles(p2d_res_zr_flt, "res", "flt");
  plotProfiles(p2d_res_zr_smt, "res", "smt");
  plotProfiles(p2d_pull_zr_prt, "pull", "prt");
  plotProfiles(p2d_pull_zr_flt, "pull", "flt");
  plotProfiles(p2d_pull_zr_smt, "pull", "smt");

  // Plotting residual
  TCanvas* residuals =
      new TCanvas("residuals", "Residual Distributions", 1200, 800);
  residuals->Divide(3, 2);

  TCanvas* pulls = new TCanvas("pulls", "Pull distributions", 1200, 800);
  pulls->Divide(3, 2);

  for (auto [vol, layers] : volLayIds) {
    for (auto lay : layers) {
      auto vlID = volLayIdCut(vol, lay)[0];

      // Residual plotting
      for (size_t ipar = 0; ipar < paramNames.size(); ipar++) {
        residuals->cd(ipar + 1);
        // Draw them
        res_smt[vlID + paramNames.at(ipar)]->Draw("");
        res_prt[vlID + paramNames.at(ipar)]->Draw("same");
        res_flt[vlID + paramNames.at(ipar)]->Draw("same");
        // Write them
        pull_smt[vlID + paramNames.at(ipar)]->Write();
        pull_prt[vlID + paramNames.at(ipar)]->Write();
        pull_flt[vlID + paramNames.at(ipar)]->Write();

        int binmax = res_smt[vlID + paramNames.at(ipar)]->GetMaximumBin();
        int bincontent =
            res_smt[vlID + paramNames.at(ipar)]->GetBinContent(binmax);

        res_smt[vlID + paramNames.at(ipar)]->GetYaxis()->SetRangeUser(
            0, bincontent * 1.2);
        TLegend* legend = new TLegend(0.7, 0.7, 0.9, 0.9);
        legend->AddEntry(res_prt[vlID + paramNames.at(ipar)], "prediction",
                         "lp");
        legend->AddEntry(res_flt[vlID + paramNames.at(ipar)], "filtering",
                         "lp");
        legend->AddEntry(res_smt[vlID + paramNames.at(ipar)], "smoothing",
                         "lp");
        legend->SetBorderSize(0);
        legend->SetFillStyle(0);
        legend->SetTextFont(42);
        legend->Draw();
      }
      if (not saveAs.empty()) {
        residuals->SaveAs((vlID + std::string("residuals.") + saveAs).c_str());
      }

      // Pull plotting & writing
      for (size_t ipar = 0; ipar < paramNames.size(); ipar++) {
        pulls->cd(ipar + 1);
        // Draw them
        pull_smt[vlID + paramNames.at(ipar)]->Draw("");
        pull_prt[vlID + paramNames.at(ipar)]->Draw("same");
        pull_flt[vlID + paramNames.at(ipar)]->Draw("same");
        // Write them
        pull_smt[vlID + paramNames.at(ipar)]->Write();
        pull_prt[vlID + paramNames.at(ipar)]->Write();
        pull_flt[vlID + paramNames.at(ipar)]->Write();

        int binmax = pull_smt[vlID + paramNames.at(ipar)]->GetMaximumBin();
        int bincontent =
            pull_smt[vlID + paramNames.at(ipar)]->GetBinContent(binmax);

        pull_smt[vlID + paramNames.at(ipar)]->GetYaxis()->SetRangeUser(
            0, bincontent * 1.2);
        TLegend* legend = new TLegend(0.7, 0.7, 0.9, 0.9);
        legend->AddEntry(pull_prt[vlID + paramNames.at(ipar)], "prediction",
                         "lp");
        legend->AddEntry(pull_flt[vlID + paramNames.at(ipar)], "filtering",
                         "lp");
        legend->AddEntry(pull_smt[vlID + paramNames.at(ipar)], "smoothing",
                         "lp");
        legend->SetBorderSize(0);
        legend->SetFillStyle(0);
        legend->SetTextFont(42);
        legend->Draw();
      }

      // Save the Canvae as pictures
      if (not saveAs.empty()) {
        pulls->SaveAs((vlID + std::string("pulls.") + saveAs).c_str());
      }
    }
  }

  // out->Close();
}


